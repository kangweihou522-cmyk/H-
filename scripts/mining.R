# Per-sample pathway activity and optional WGCNA; no small-sample network inference.
trait_matrix <- function(m,selected) {
 selected<-unique(c('group',unlist(selected)))
 columns<-list()
 for(name in selected){
  x<-m[[name]]
  if(is.numeric(x) && length(unique(x))>1)columns[[name]]<-x
  else if(length(unique(x))>1){
   for(level in sort(unique(as.character(x))))columns[[paste(name,level,sep='=')]]<-as.numeric(as.character(x)==level)
  }
 }
 if(!length(columns))return(matrix(numeric(),nrow=nrow(m),ncol=0))
 out<-as.matrix(as.data.frame(columns,check.names=FALSE));rownames(out)<-m$sample;out
}
contrast_vector <- function(m,design,cmp) {
 template<-m[c(1,1),,drop=FALSE]
 template$group<-factor(c(cmp$treatment,cmp$control),levels=levels(m$group))
 x<-model.matrix(design,template);x[1,]-x[2,]
}
run_mining <- function(vs,norm,m,mm,annotation,root,cfg) {
 mining<-list(pathway_activity=list(state='skipped',reason='No GMT supplied or activity disabled'),wgcna=list(state='skipped',reason='Not selected'))
 traits<-trait_matrix(m,cfg$traits)
 write.csv(traits,file.path(root,'traits.csv'))
 if(cfg$pathway_activity && length(annotation$gmt)){
  ids<-if(cfg$gmt_id=='ENSEMBL')annotation$mapped$id else if(cfg$gmt_id=='SYMBOL')annotation$mapped$symbol else annotation$mapped$entrez
  row_order<-order(rowMeans(norm),decreasing=TRUE)
  keep<-row_order[!is.na(ids[row_order]) & nzchar(ids[row_order]) & !duplicated(ids[row_order])]
  expr<-vs[keep,,drop=FALSE];rownames(expr)<-ids[keep]
  activity<-GSVA::gsva(expr,annotation$gmt,method='ssgsea',min.sz=cfg$min_set,max.sz=cfg$max_set,ssgsea.norm=TRUE,parallel.sz=1,verbose=FALSE)
  if(nrow(activity)){
   write.csv(activity,file.path(root,'pathway-activity.csv'))
   fit<-limma::lmFit(activity,mm)
   design<-reformulate(c(vapply(cfg$covariates,`[[`,character(1),'name'),'group'))
   for(cmp in cfg$comparisons){
    vector<-contrast_vector(m,design,cmp)
    model<-limma::eBayes(limma::contrasts.fit(fit,matrix(vector,ncol=1)),trend=TRUE)
    r<-limma::topTable(model,coef=1,number=Inf,sort.by='P');r$term<-rownames(r)
    r$description<-unname(annotation$gmt_desc[r$term]);write.csv(r,file.path(root,'contrasts',cmp$name,'pathway-activity-differential.csv'),row.names=FALSE)
   }
   if(ncol(traits)){
    cor<-cor(t(activity),traits);p<-WGCNA::corPvalueStudent(cor,nrow(m))
    rows<-expand.grid(term=rownames(cor),trait=colnames(cor),stringsAsFactors=FALSE);rows$correlation<-as.vector(cor);rows$pvalue<-as.vector(p);rows$padj<-p.adjust(rows$pvalue,'BH')
    write.csv(rows,file.path(root,'pathway-trait.csv'),row.names=FALSE)
   }
   mining$pathway_activity<-list(state='completed',method='GSVA ssGSEA; VST expression; limma moderated test with the same covariates',sets=nrow(activity))
  }else mining$pathway_activity<-list(state='skipped',reason='No gene sets passed size/mapping filters')
 }
 if(cfg$wgcna){
  independent<-if('subject'%in%names(m))length(unique(m$subject)) else nrow(m)
  if(independent<nrow(m)){
   mining$wgcna<-list(state='skipped',reason='Repeated subjects detected. WGCNA currently requires one sample per independent subject; repeated measures need separate modelling.')
  }else if(independent<15){
   mining$wgcna<-list(state='skipped',reason=paste('WGCNA requires at least 15 independent samples; 20 or more recommended. Available:',independent))
  }else{
   var<-apply(vs,1,var);choose<-head(order(var,decreasing=TRUE),cfg$wgcna_genes);choose<-choose[var[choose]>0]
   dat<-t(vs[choose,,drop=FALSE]);check<-WGCNA::goodSamplesGenes(dat,verbose=0)
   if(!check$allOK)stop('WGCNA goodSamplesGenes failed; inspect input before excluding samples')
   suppressPackageStartupMessages(library(WGCNA))
   if(cfg$threads>=2)WGCNA::allowWGCNAThreads(nThreads=cfg$threads) else WGCNA::disableWGCNAThreads()
   soft<-WGCNA::pickSoftThreshold(dat,powerVector=c(1:10,12,14,16,18,20),networkType='signed',verbose=0)
   write.csv(soft$fitIndices,file.path(root,'wgcna-soft-threshold.csv'),row.names=FALSE)
   net<-WGCNA::blockwiseModules(dat,power=cfg$wgcna_power,networkType='signed',TOMType='signed',minModuleSize=30,maxBlockSize=3000,mergeCutHeight=.25,nThreads=cfg$threads,verbose=0,randomSeed=42)
   eig<-WGCNA::orderMEs(net$MEs);write.csv(eig,file.path(root,'wgcna-eigengenes.csv'))
   kme<-WGCNA::signedKME(dat,eig)
   membership<-vapply(seq_len(ncol(dat)),function(i){col<-match(paste0('kME',net$colors[i]),colnames(kme));if(is.na(col))NA_real_ else kme[i,col]},numeric(1))
   modules<-data.frame(gene=colnames(dat),module=net$colors,kME=membership);write.csv(modules,file.path(root,'wgcna-modules.csv'),row.names=FALSE)
   if(ncol(traits)){
    cor<-cor(eig,traits);p<-WGCNA::corPvalueStudent(cor,nrow(m));rows<-expand.grid(module=rownames(cor),trait=colnames(cor),stringsAsFactors=FALSE);rows$correlation<-as.vector(cor);rows$pvalue<-as.vector(p);rows$padj<-p.adjust(rows$pvalue,'BH');write.csv(rows,file.path(root,'wgcna-module-trait.csv'),row.names=FALSE)
   }
   hubs<-head(modules$gene[order(abs(modules$kME),decreasing=TRUE)],120);edges<-cor(dat[,hubs,drop=FALSE]);adj<-((1+edges)/2)^cfg$wgcna_power;idx<-which(upper.tri(adj)&adj>.2,arr.ind=TRUE)
   if(nrow(idx)){out<-data.frame(source=rownames(adj)[idx[,1]],target=colnames(adj)[idx[,2]],weight=adj[idx]);out<-head(out[order(out$weight,decreasing=TRUE),],500)}else out<-data.frame(source=character(),target=character(),weight=numeric())
   write.csv(out,file.path(root,'wgcna-network.csv'),row.names=FALSE)
   mining$wgcna<-list(state='completed',method='WGCNA signed network; VST; top variable genes; raw module/trait Pearson associations; BH',samples=nrow(dat),genes=ncol(dat),power=cfg$wgcna_power,modules=length(unique(net$colors)))
  }
 }
 mining$regulatory_activity<-regulatory_activity(vs,norm,m,mm,annotation,root,cfg)
 write_json(mining,file.path(root,'mining-status.json'),auto_unbox=TRUE,pretty=TRUE)
}

regulatory_activity <- function(vs,norm,m,mm,annotation,root,cfg){
 path<-file.path(root,'regulatory-network.tsv')
 if(!file.exists(path))return(list(state='skipped',reason='Need a mouse source/target/weight regulatory network, e.g. an authorized DoRothEA mouse export.'))
 ids<-if(cfg$tf_id=='ENSEMBL')annotation$mapped$id else if(cfg$tf_id=='SYMBOL')annotation$mapped$symbol else annotation$mapped$entrez
 order<-order(rowMeans(norm),decreasing=TRUE);keep<-order[!is.na(ids[order]) & nzchar(ids[order]) & !duplicated(ids[order])]
 expr<-vs[keep,,drop=FALSE];rownames(expr)<-ids[keep]
 network<-read.delim(path,stringsAsFactors=FALSE)
 mapped<-network[network$target%in%rownames(expr),]
 eligible<-names(which(table(mapped$source)>=5));mapped<-mapped[mapped$source%in%eligible,]
 if(!nrow(mapped))return(list(state='skipped',reason='No regulator has at least 5 detected target genes. Check mouse identifiers and network.'))
 estimate<-decoupleR::run_ulm(expr,mapped,.source='source',.target='target',.mor='weight',minsize=5)
 scores<-matrix(NA_real_,nrow=length(unique(estimate$source)),ncol=nrow(m),dimnames=list(unique(estimate$source),m$sample));scores[cbind(match(estimate$source,rownames(scores)),match(estimate$condition,colnames(scores)))]<-estimate$score
 write.csv(scores,file.path(root,'tf-activity.csv'));write.csv(estimate,file.path(root,'tf-activity-details.csv'),row.names=FALSE)
 fit<-limma::lmFit(scores,mm);design<-reformulate(c(vapply(cfg$covariates,`[[`,character(1),'name'),'group'))
 for(cmp in cfg$comparisons){
  vector<-contrast_vector(m,design,cmp);model<-limma::eBayes(limma::contrasts.fit(fit,matrix(vector,ncol=1)),trend=TRUE)
  out<-limma::topTable(model,coef=1,number=Inf,sort.by='P');out$term<-rownames(out);out$description<-out$term
  write.csv(out,file.path(root,'contrasts',cmp$name,'tf-activity-differential.csv'),row.names=FALSE)
 }
 list(state='completed',method='decoupleR ULM; mouse prior network; VST; limma contrast with the same covariates',regulators=nrow(scores),targets=length(unique(mapped$target)))
}
