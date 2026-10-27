# Standalone analysis from complete DE results or a declared transformed matrix.
suppressPackageStartupMessages({library(jsonlite);library(GO.db);library(fgsea)})
a<-commandArgs(TRUE);root<-a[1];cfg<-fromJSON(file.path(root,'config.json'),simplifyVector=FALSE)
script<-sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE))
script_dir<-Sys.getenv('RNASEQ_SCRIPT_DIR',unset=dirname(gsub('~+~',' ',script,fixed=TRUE)));source(file.path(script_dir,'enrichment.R'));source(file.path(script_dir,'mining.R'))
stage<-cfg$analysis_stage;m<-read.csv(file.path(root,'samples.csv'),colClasses='character',check.names=FALSE)
for(cmp in cfg$comparisons)dir.create(file.path(root,'contrasts',cmp$name),recursive=TRUE,showWarnings=FALSE)
summary<-list(design='No DESeq2 fit in this standalone module',samples=nrow(m),contrasts=list(),analysis_stage=stage,expression_unit=cfg$expression_unit)
vs<-NULL
if(file.exists(file.path(root,'expression-input.csv'))){
 rownames(m)<-m$sample;m$group<-factor(m$group)
 if(cfg$control %in% levels(m$group))m$group<-relevel(m$group,cfg$control)
 for(c in cfg$covariates)m[[c$name]]<-if(c$type=='numeric')as.numeric(m[[c$name]]) else factor(m[[c$name]])
 vs<-as.matrix(read.csv(file.path(root,'expression-input.csv'),row.names=1,check.names=FALSE,colClasses=c('character',rep('numeric',nrow(m)))))
 vs<-vs[,m$sample,drop=FALSE]
 # "vst.csv" is the common plotting interface; the declared input unit is saved
 # and displayed. No VST, size-factor estimation or fake count data is generated.
 write.csv(vs,file.path(root,'vst.csv'))
 variables<-c(vapply(cfg$covariates,`[[`,character(1),'name'),if(nlevels(m$group)>1)'group')
 design<-if(length(variables))reformulate(variables) else ~1
 mm<-model.matrix(design,m);write.csv(mm,file.path(root,'design-matrix.csv'))
 genes<-head(order(apply(vs,1,var),decreasing=TRUE),500)
 pc<-prcomp(t(vs[genes,,drop=FALSE]),center=TRUE,scale.=FALSE)
 write.csv(data.frame(sample=m$sample,group=m$group,pc$x[,seq_len(min(5,ncol(pc$x))),drop=FALSE]),file.path(root,'pca-scores.csv'),row.names=FALSE)
 write.csv(data.frame(PC=seq_along(pc$sdev),variance=pc$sdev^2/sum(pc$sdev^2)),file.path(root,'pca-variance.csv'),row.names=FALSE)
 write.csv(cor(vs),file.path(root,'sample-correlation.csv'));write.csv(as.matrix(dist(t(vs))),file.path(root,'sample-distance.csv'))
 summary$design<-paste(deparse(design),collapse='');summary$pca_genes<-length(genes);summary$genes_input<-nrow(vs);summary$genes_retained<-nrow(vs)
}
if(stage=='enrichment'){
 r<-read.csv(file.path(root,'differential-input.csv'),colClasses=c(gene='character'),check.names=FALSE)
 summary$padj_source<-if('padj'%in%names(r))'Imported padj unchanged' else 'BH across all supplied finite pvalue (complete tested table required)'
 if(!'padj'%in%names(r)){r$padj<-NA_real_;valid<-is.finite(r$pvalue);r$padj[valid]<-p.adjust(r$pvalue[valid],'BH')}
 if(!'baseMean'%in%names(r))r$baseMean<-NA_real_
 annotation<-load_annotation(r$gene,root,cfg)
 if(!is.null(annotation$warning))writeLines(annotation$warning,file.path(root,'annotation-warning.txt'))
 write_json(list(go=annotation$go,pathway=annotation$gmt,go_desc=as.list(annotation$go_desc),pathway_desc=as.list(annotation$gmt_desc)),file.path(root,'gene-sets.json'),auto_unbox=TRUE)
 if(!'symbol'%in%names(r))r$symbol<-annotation$mapped$symbol[match(r$gene,annotation$mapped$gene)]
 cmp<-cfg$comparisons[[1]];directory<-file.path(root,'contrasts',cmp$name);dir.create(directory,recursive=TRUE,showWarnings=FALSE)
 counts<-list()
 for(profile in c('nominal','strict')){
  column<-if(profile=='nominal')'pvalue' else 'padj';lfc<-if(profile=='nominal')0 else 1
  rp<-r;rp$significant<-is.finite(rp[[column]]) & rp[[column]]<.05 & is.finite(rp$log2FoldChange) & abs(rp$log2FoldChange)>lfc
  rp$direction<-ifelse(rp$significant,ifelse(rp$log2FoldChange>0,'up','down'),'not_significant')
  pp<-file.path(directory,'profiles',profile);dir.create(pp,recursive=TRUE,showWarnings=FALSE)
  write.csv(rp,file.path(pp,'differential.csv'),row.names=FALSE)
  for(direction in c('all','up','down'))write.csv(rp[rp$significant & (direction=='all'|rp$direction==direction),],file.path(pp,paste0('deg-',direction,'.csv')),row.names=FALSE)
  run_enrichment(rp,annotation,pp,cfg,run_gsea=profile=='strict' && cfg$run_gsea)
  counts[[profile]]<-list(profile=profile,column=column,alpha=.05,lfc=lfc,total=sum(rp$significant),up=sum(rp$direction=='up'),down=sum(rp$direction=='down'))
 }
 active<-counts[[cfg$deg_profile]]
 for(file in c('differential.csv','deg-all.csv','deg-up.csv','deg-down.csv','go.csv','pathway.csv'))file.copy(file.path(directory,'profiles',cfg$deg_profile,file),file.path(directory,file),overwrite=TRUE)
 for(file in c('gsea.csv','gsea-ranks.csv','pathway_gsea.csv','pathway_gsea-ranks.csv'))if(file.exists(file.path(directory,'profiles','strict',file)))file.copy(file.path(directory,'profiles','strict',file),file.path(directory,file),overwrite=TRUE)
 write_json(active,file.path(directory,'active-profile.json'),auto_unbox=TRUE)
 summary$contrasts[[cmp$name]]<-list(name=cmp$name,control=cmp$control,treatment=cmp$treatment,tested=nrow(r),up=active$up,down=active$down,significance=active$column,alpha=.05,lfc=active$lfc,active_profile=cfg$deg_profile,profiles=counts,shrinkage='Imported LFC unchanged; no DESeq2 refit')
 summary$genes_retained<-nrow(r);summary$gsea_engine<-if(cfg$run_gsea)'fgsea multilevel; imported complete stat ranking' else 'Skipped: no requested / sufficient stat ranking'
 for(file in list.files(directory,full.names=TRUE)[!dir.exists(list.files(directory,full.names=TRUE))])file.copy(file,file.path(root,basename(file)),overwrite=TRUE)
 mining<-list(pathway_activity=list(state='skipped',reason='Standalone enrichment module'),wgcna=list(state='skipped',reason='Standalone enrichment module'))
 write_json(mining,file.path(root,'mining-status.json'),auto_unbox=TRUE,pretty=TRUE)
}else{
 for(cmp in cfg$comparisons)summary$contrasts[[cmp$name]]<-list(name=cmp$name,control=cmp$control,treatment=cmp$treatment,method='Activity score comparison; no gene DE fit')
 annotation<-load_annotation(rownames(vs),root,cfg)
 if(!is.null(annotation$warning))writeLines(annotation$warning,file.path(root,'annotation-warning.txt'))
 run_mining(vs,vs,m,mm,annotation,root,cfg)
 if(!is.null(cfg$expression_unit) && cfg$expression_unit!='vst'){
  mining<-read_json(file.path(root,'mining-status.json'),simplifyVector=FALSE)
  for(key in names(mining))if(!is.null(mining[[key]]$method))mining[[key]]$method<-gsub('VST',cfg$expression_unit,mining[[key]]$method,fixed=TRUE)
  write_json(mining,file.path(root,'mining-status.json'),auto_unbox=TRUE,pretty=TRUE)
 }
}
write_json(summary,file.path(root,'analysis-summary.json'),auto_unbox=TRUE,pretty=TRUE)
writeLines(capture.output(sessionInfo()),file.path(root,'R-session.txt'))
