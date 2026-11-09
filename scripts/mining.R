# Per-sample pathway activity and optional WGCNA; no small-sample network inference.
trait_matrix <- function(m,selected) {
 selected<-unique(c('group',unlist(selected)))
 columns<-list()
 for(name in selected){
  x<-m[[name]]
  numeric<-suppressWarnings(as.numeric(as.character(x)))
  if(name!='group' && all(is.finite(numeric)))x<-numeric
  if(is.numeric(x) && length(unique(x))>1)columns[[name]]<-x
  else if(length(unique(x))>1){
   for(level in sort(unique(as.character(x))))columns[[paste(name,level,sep='=')]]<-as.numeric(as.character(x)==level)
  }
 }
 if(!length(columns))return(matrix(numeric(),nrow=nrow(m),ncol=0))
 out<-as.matrix(as.data.frame(columns,check.names=FALSE));rownames(out)<-m$sample;out
}
select_soft_power <- function(fits,mode,manual) {
 if(is.null(mode)||mode=='manual')return(manual)
 pass<-fits$Power[is.finite(fits$SFT.R.sq) & is.finite(fits$slope) & fits$slope<0 & fits$SFT.R.sq>=.85]
 if(length(pass))min(pass) else NA_real_
}
module_enrichment <- function(modules,norm,annotation,root,cfg){
 # Universe is all genes entered into this network (including unassigned grey),
 # not the whole genome or only DEGs. Adjust across module x term within ontology.
 modules<-modules[order(-rowMeans(norm[modules$gene,,drop=FALSE]),modules$gene),]
 map<-annotation$mapped[match(modules$gene,annotation$mapped$gene),]
 outputs<-list()
 for(database in c('go','pathway')){
  ids<-if(database=='go'||cfg$gmt_id=='ENSEMBL')map$id else if(cfg$gmt_id=='SYMBOL')map$symbol else map$entrez
  valid<-!is.na(ids)&nzchar(ids);valid[valid]<-!duplicated(ids[valid])
  sets<-if(database=='go')annotation$go else annotation$gmt
  desc<-if(database=='go')annotation$go_desc else annotation$gmt_desc
  ont<-if(database=='go')annotation$ontology else character()
  for(module in setdiff(unique(modules$module),'grey')){
   x<-ora(sets,desc,ids[valid],ids[valid & modules$module==module],ont,cfg)
   if(nrow(x)){x$module<-module;x$database<-database;outputs[[length(outputs)+1]]<-x}
  }
 }
 if(length(outputs)){
  out<-do.call(rbind,outputs);out$padj<-ave(out$pvalue,interaction(out$database,out$ontology,drop=TRUE),FUN=function(x)p.adjust(x,'BH'))
  write.csv(out,file.path(root,'wgcna-module-enrichment.csv'),row.names=FALSE)
 }
}
run_mining_impl <- function(vs,norm,m,mm,annotation,root,cfg) {
 mining<-list(pathway_activity=list(state='skipped',reason='No GMT supplied or activity disabled'),wgcna=list(state='skipped',reason='Not selected'))
 traits<-trait_matrix(m,cfg$traits)
 independent<-if(cfg$subject_column%in%names(m))length(unique(m[[cfg$subject_column]])) else nrow(m)
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
   design<-model_formula(cfg)
   for(cmp in cfg$comparisons){
    vector<-contrast_vector(m,design,cmp,cfg)
    model<-limma::eBayes(limma::contrasts.fit(fit,matrix(vector,ncol=1)),trend=TRUE)
    r<-limma::topTable(model,coef=1,number=Inf,sort.by='P');r$term<-rownames(r)
    r$description<-unname(annotation$gmt_desc[r$term]);write.csv(r,file.path(root,'contrasts',cmp$name,'pathway-activity-differential.csv'),row.names=FALSE)
   }
   if(ncol(traits) && independent==nrow(m)){
    rows<-feature_traits(activity,traits,m,mm,cfg,'term')
    if(nrow(rows))write.csv(rows,file.path(root,'pathway-trait.csv'),row.names=FALSE)
   }
   mining$pathway_activity<-list(state='completed',method='GSVA ssGSEA; VST expression; limma moderated test with the same covariates',sets=nrow(activity))
   if(independent<nrow(m))mining$pathway_activity$trait_association<-'Skipped raw trait association: repeated subjects require a repeated-measures model.'
  }else mining$pathway_activity<-list(state='skipped',reason='No gene sets passed size/mapping filters')
 }
 if(cfg$wgcna){
  independent<-if(cfg$subject_column%in%names(m))length(unique(m[[cfg$subject_column]])) else nrow(m)
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
   power<-select_soft_power(soft$fitIndices,cfg$wgcna_power_mode,cfg$wgcna_power)
   if(!is.finite(power)){
     mining$wgcna<-list(state='skipped',reason='No power reached signed scale-free fit R2 >= 0.85 with negative slope. Review diagnostics; choose a justified manual power in a new job.')
     mining$regulatory_activity<-regulatory_activity(vs,norm,m,mm,annotation,root,cfg)
     write_json(mining,file.path(root,'mining-status.json'),auto_unbox=TRUE,pretty=TRUE);return(invisible(NULL))
   }
   net<-WGCNA::blockwiseModules(dat,power=power,networkType='signed',TOMType='signed',corType=if(is.null(cfg$wgcna_correlation))'pearson' else cfg$wgcna_correlation,maxPOutliers=.05,minModuleSize=30,maxBlockSize=3000,mergeCutHeight=.25,nThreads=cfg$threads,verbose=0,randomSeed=42)
   if(all(net$colors=='grey')){
    write.csv(data.frame(gene=colnames(dat),module='grey',kME=NA_real_),file.path(root,'wgcna-modules.csv'),row.names=FALSE)
    mining$wgcna<-list(state='skipped',reason='No non-grey coexpression modules detected; diagnostics retained.')
    mining$regulatory_activity<-regulatory_activity(vs,norm,m,mm,annotation,root,cfg)
    write_json(mining,file.path(root,'mining-status.json'),auto_unbox=TRUE,pretty=TRUE);return(invisible(NULL))
   }
   eig<-WGCNA::orderMEs(net$MEs);eig<-eig[,colnames(eig)!='MEgrey',drop=FALSE];write.csv(eig,file.path(root,'wgcna-eigengenes.csv'))
   kme<-WGCNA::signedKME(dat,eig)
   membership<-vapply(seq_len(ncol(dat)),function(i){col<-match(paste0('kME',net$colors[i]),colnames(kme));if(is.na(col)||net$colors[i]=='grey')NA_real_ else kme[i,col]},numeric(1))
   modules<-data.frame(gene=colnames(dat),module=net$colors,kME=membership);write.csv(modules,file.path(root,'wgcna-modules.csv'),row.names=FALSE)
   module_enrichment(modules,norm,annotation,root,cfg)
   network_stability(dat,modules,root,cfg)
   if(ncol(traits)){
    rows<-feature_traits(t(eig),traits,m,mm,cfg,'module');if(nrow(rows))write.csv(rows,file.path(root,'wgcna-module-trait.csv'),row.names=FALSE)
   }
   hubs<-head(modules$gene[order(abs(modules$kME),decreasing=TRUE,na.last=NA)],120);edges<-cor(dat[,hubs,drop=FALSE]);adj<-((1+edges)/2)^power;idx<-which(upper.tri(adj)&adj>.2,arr.ind=TRUE)
   if(nrow(idx)){out<-data.frame(source=rownames(adj)[idx[,1]],target=colnames(adj)[idx[,2]],weight=adj[idx]);out<-head(out[order(out$weight,decreasing=TRUE),],500)}else out<-data.frame(source=character(),target=character(),weight=numeric())
   write.csv(out,file.path(root,'wgcna-network.csv'),row.names=FALSE)
   mining$wgcna<-list(state='completed',method='WGCNA signed network; VST; top variable genes; raw and covariate-adjusted module/trait associations; BH within method',samples=nrow(dat),genes=ncol(dat),power=power,power_mode=if(is.null(cfg$wgcna_power_mode))'manual' else cfg$wgcna_power_mode,modules=length(setdiff(unique(net$colors),'grey')))
  }
 }
 mining$regulatory_activity<-regulatory_activity(vs,norm,m,mm,annotation,root,cfg)
 write_json(mining,file.path(root,'mining-status.json'),auto_unbox=TRUE,pretty=TRUE)
}

regulatory_activity <- function(vs,norm,m,mm,annotation,root,cfg){
 if(identical(cfg$tf_enabled,FALSE))return(list(state='skipped',reason='Not selected in this module'))
 prefix<-if(is.null(cfg$activity_prefix))'tf' else cfg$activity_prefix
 path<-file.path(root,if(prefix=='tf')'regulatory-network.tsv' else 'pathway-network.tsv')
 if(!file.exists(path))return(list(state='skipped',reason=if(prefix=='tf')'Need a mouse source/target/weight regulatory network, e.g. an authorized DoRothEA mouse export.' else 'Need a real mouse response footprint model; no built-in PROGENy database.'))
 ids<-if(cfg$tf_id=='ENSEMBL')annotation$mapped$id else if(cfg$tf_id=='SYMBOL')annotation$mapped$symbol else annotation$mapped$entrez
 order<-order(rowMeans(norm),decreasing=TRUE);keep<-order[!is.na(ids[order]) & nzchar(ids[order]) & !duplicated(ids[order])]
 expr<-vs[keep,,drop=FALSE];rownames(expr)<-ids[keep]
 network<-read.delim(path,stringsAsFactors=FALSE)
 mapped<-network[network$target%in%rownames(expr),]
 eligible<-names(which(table(mapped$source)>=5));mapped<-mapped[mapped$source%in%eligible,]
 if(!nrow(mapped))return(list(state='skipped',reason='No regulator has at least 5 detected target genes. Check mouse identifiers and network.'))
 estimate<-decoupleR::run_ulm(expr,mapped,.source='source',.target='target',.mor='weight',minsize=5)
 scores<-matrix(NA_real_,nrow=length(unique(estimate$source)),ncol=nrow(m),dimnames=list(unique(estimate$source),m$sample));scores[cbind(match(estimate$source,rownames(scores)),match(estimate$condition,colnames(scores)))]<-estimate$score
 write.csv(scores,file.path(root,paste0(prefix,'-activity.csv')));write.csv(estimate,file.path(root,paste0(prefix,'-activity-details.csv')),row.names=FALSE)
 fit<-limma::lmFit(scores,mm);design<-model_formula(cfg)
 for(cmp in cfg$comparisons){
  vector<-contrast_vector(m,design,cmp,cfg);model<-limma::eBayes(limma::contrasts.fit(fit,matrix(vector,ncol=1)),trend=TRUE)
  out<-limma::topTable(model,coef=1,number=Inf,sort.by='P');out$term<-rownames(out);out$description<-out$term
  write.csv(out,file.path(root,'contrasts',cmp$name,paste0(prefix,'-activity-differential.csv')),row.names=FALSE)
 }
 list(state='completed',method='decoupleR ULM; mouse prior network; VST; limma contrast with the same covariates',regulators=nrow(scores),targets=length(unique(mapped$target)),model_source=if(prefix=='tf')cfg$tf_network else cfg$pathway_source)
}

# Optional analyses have separate failure states; primary DE survives a mining error.
run_mining <- function(vs,norm,m,mm,annotation,root,cfg){
 if(is.null(cfg$subject_column))cfg$subject_column<-'subject'
 state<-list()
 for(key in c('pathway_activity','wgcna','regulatory_activity','footprint_activity')){
  sub<-cfg;sub$pathway_activity<-key=='pathway_activity'&&isTRUE(cfg$pathway_activity);sub$wgcna<-key=='wgcna'&&isTRUE(cfg$wgcna);sub$tf_enabled<-key%in%c('regulatory_activity','footprint_activity')
  if(key=='footprint_activity'){sub$activity_prefix<-'footprint';sub$tf_id<-if(is.null(cfg$pathway_network_id))'ENSEMBL' else cfg$pathway_network_id}
  state[[key]]<-tryCatch({run_mining_impl(vs,norm,m,mm,annotation,root,sub);fromJSON(file.path(root,'mining-status.json'),simplifyVector=FALSE)[[if(key=='footprint_activity')'regulatory_activity' else key]]},error=function(e)list(state='failed',reason=conditionMessage(e)))
  atomic_json(state,file.path(root,'mining-status.json'))
 }
 invisible(state)
}
network_stability <- function(dat,modules,root,cfg){
 B<-if(is.null(cfg$wgcna_resamples))0 else cfg$wgcna_resamples
 if(B<1)return(invisible(NULL))
 # Sample subsampling checks membership stability of fixed modules, not re-discovery or external preservation.
 if(nrow(dat)<20){atomic_json(list(state='skipped',reason='Stability needs at least20 independent samples'),file.path(root,'wgcna-stability-status.json'));return(invisible(NULL))}
 set.seed(42);out<-list()
 for(b in seq_len(B)){
  idx<-sample(seq_len(nrow(dat)),max(15,floor(.8*nrow(dat))))
  for(module in setdiff(unique(modules$module),'grey')){
   genes<-modules$gene[modules$module==module];x<-dat[idx,genes,drop=FALSE];pc<-prcomp(x,center=TRUE,scale.=FALSE)$x[,1]
   km<-abs(as.vector(cor(x,pc)));original<-abs(modules$kME[match(genes,modules$gene)])
   out[[length(out)+1]]<-data.frame(module=module,resample=b,samples=length(idx),membership_rank_rho=suppressWarnings(cor(km,original,method='spearman')),top_hub_jaccard={a<-head(genes[order(km,decreasing=TRUE)],10);z<-head(genes[order(original,decreasing=TRUE)],10);length(intersect(a,z))/length(union(a,z))})
  }
 }
 if(length(out))write.csv(do.call(rbind,out),file.path(root,'wgcna-stability.csv'),row.names=FALSE)
 atomic_json(list(state='completed',method='80% sample subsampling; fixed-module membership rank and top10 hub stability; not module preservation',resamples=B),file.path(root,'wgcna-stability-status.json'))
}
