suppressPackageStartupMessages({library(DESeq2);library(tximport);library(GO.db);library(jsonlite);library(fgsea)})
a <- commandArgs(TRUE);root<-a[1];cfg<-fromJSON(file.path(root,'config.json'),simplifyVector=FALSE)
if(is.null(cfg$deg_profile))cfg$deg_profile<-'strict'
script <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE))
# Rscript encodes spaces as ~+~ in --file. The app supplies the exact mounted path.
script_dir<-Sys.getenv('RNASEQ_SCRIPT_DIR',unset=dirname(gsub('~+~',' ',script,fixed=TRUE)))
source(file.path(script_dir,'enrichment.R'));source(file.path(script_dir,'mining.R'))
m<-read.csv(file.path(root,'samples.csv'),check.names=FALSE,colClasses='character');rownames(m)<-m$sample
m$group<-relevel(factor(m$group),cfg$control)
for(c in cfg$covariates)m[[c$name]]<-if(c$type=='numeric')as.numeric(m[[c$name]]) else factor(m[[c$name]])
design<-reformulate(c(vapply(cfg$covariates,`[[`,character(1),'name'),'group'))
mm<-model.matrix(design,m);if(qr(mm)$rank<ncol(mm)||nrow(mm)<=ncol(mm))stop('Design is confounded or lacks residual degrees of freedom.')
if(file.exists(file.path(root,'counts.csv'))){
 counts<-read.csv(file.path(root,'counts.csv'),row.names=1,check.names=FALSE,colClasses=c('character',rep('numeric',nrow(m))))
 dds<-DESeqDataSetFromMatrix(as.matrix(counts[,m$sample,drop=FALSE]),m,design)
}else{
 files<-file.path(root,'quant',m$sample,'quant.sf');names(files)<-m$sample
 map<-read.delim(file.path(root,'tx2gene.tsv'),header=TRUE,colClasses='character')
 map[[1]]<-sub('\\..*$','',sub('\\|.*$','',map[[1]]))
 txi<-tximport(files,type='salmon',tx2gene=map[,1:2],ignoreTxVersion=TRUE,ignoreAfterBar=TRUE)
 write.csv(txi$abundance,file.path(root,'tpm.csv'))
 dds<-DESeqDataSetFromTximport(txi,m,design)
}
raw_all<-counts(dds);before<-nrow(dds);dds<-dds[rowSums(counts(dds)>=cfg$min_count)>=cfg$min_samples,];if(nrow(dds)<10)stop('Too few expressed genes after filtering')
dds<-DESeq(dds,quiet=TRUE)
write.csv(mm,file.path(root,'design-matrix.csv'))
norm<-counts(dds,normalized=TRUE);vs<-assay(varianceStabilizingTransformation(dds,blind=FALSE))
write.csv(raw_all,file.path(root,'raw-counts.csv'));write.csv(norm,file.path(root,'normalized.csv'));write.csv(vs,file.path(root,'vst.csv'))
pca_genes<-order(apply(vs,1,var),decreasing=TRUE)[seq_len(min(500,nrow(vs)))]
pca<-prcomp(t(vs[pca_genes,,drop=FALSE]),center=TRUE,scale.=FALSE)
write.csv(data.frame(sample=m$sample,group=m$group,pca$x[,seq_len(min(5,ncol(pca$x))),drop=FALSE]),file.path(root,'pca-scores.csv'),row.names=FALSE)
write.csv(data.frame(PC=seq_along(pca$sdev),variance=pca$sdev^2/sum(pca$sdev^2)),file.path(root,'pca-variance.csv'),row.names=FALSE)
write.csv(cor(vs),file.path(root,'sample-correlation.csv'));write.csv(as.matrix(dist(t(vs))),file.path(root,'sample-distance.csv'))
qc<-data.frame(sample=m$sample,group=m$group,library_size=colSums(raw_all),detected_genes=colSums(raw_all>0),median_normalized=apply(norm,2,median),size_factor=if(is.null(sizeFactors(dds)))NA_real_ else sizeFactors(dds))
write.csv(qc,file.path(root,'sample-qc.csv'),row.names=FALSE)
disp<-data.frame(gene=rownames(dds),baseMean=rowMeans(norm),dispersion=dispersions(dds),dispersion_gene=mcols(dds)$dispGeneEst,dispersion_fit=mcols(dds)$dispFit,cooks_max=apply(assays(dds)[['cooks']],1,max,na.rm=TRUE));write.csv(disp,file.path(root,'dispersion.csv'),row.names=FALSE)
annotation<-load_annotation(rownames(dds),root,cfg)
if(!is.null(annotation$warning))writeLines(annotation$warning,file.path(root,'annotation-warning.txt'))
write_json(list(go=annotation$go,pathway=annotation$gmt,go_desc=as.list(annotation$go_desc),pathway_desc=as.list(annotation$gmt_desc)),file.path(root,'gene-sets.json'),auto_unbox=TRUE)
summary<-list()
for(cmp in cfg$comparisons){
 directory<-file.path(root,'contrasts',cmp$name);dir.create(directory,recursive=TRUE,showWarnings=FALSE)
 result<-results(dds,contrast=c('group',cmp$treatment,cmp$control),alpha=cfg$alpha)
 r<-as.data.frame(result);r$gene<-rownames(r);r$log2FoldChange_raw<-r$log2FoldChange
 if(cfg$shrinkage=='normal'){
  shrunk<-lfcShrink(dds,contrast=c('group',cmp$treatment,cmp$control),res=result,type='normal',quiet=TRUE)
  r$log2FoldChange<-shrunk$log2FoldChange;r$lfcSE<-shrunk$lfcSE
 }
 r$symbol<-annotation$mapped$symbol[match(r$gene,annotation$mapped$gene)]
 profiles<-list(nominal=list(column='pvalue',alpha=.05,lfc=0),strict=list(column='padj',alpha=.05,lfc=1))
 profile_counts<-list()
 r$mean_control<-rowMeans(norm[,m$group==cmp$control,drop=FALSE]);r$mean_treatment<-rowMeans(norm[,m$group==cmp$treatment,drop=FALSE])
 r<-r[order(r$padj,na.last=TRUE),]
 # fgsea is run once on all retained ranked genes, independent of DEG profile.
 for(profile in names(profiles)){
  spec<-profiles[[profile]];rp<-r
  rp$significant<-is.finite(rp[[spec$column]]) & is.finite(rp$log2FoldChange) & rp[[spec$column]]<spec$alpha & abs(rp$log2FoldChange)>spec$lfc
  rp$direction<-ifelse(rp$significant,ifelse(rp$log2FoldChange>0,'up','down'),'not_significant')
  pp<-file.path(directory,'profiles',profile);dir.create(pp,recursive=TRUE,showWarnings=FALSE)
  write.csv(rp,file.path(pp,'differential.csv'),row.names=FALSE)
  for(direction in c('all','up','down'))write.csv(rp[rp$significant & (direction=='all' | rp$direction==direction),],file.path(pp,paste0('deg-',direction,'.csv')),row.names=FALSE)
  run_enrichment(rp,annotation,pp,cfg,run_gsea=profile=='strict')
  profile_counts[[profile]]<-list(profile=profile,column=spec$column,alpha=spec$alpha,lfc=spec$lfc,operator='strict < / >',total=sum(rp$significant),up=sum(rp$direction=='up'),down=sum(rp$direction=='down'))
 }
 active<-cfg$deg_profile
 pp<-file.path(directory,'profiles',active)
 for(file in c('differential.csv','deg-all.csv','deg-up.csv','deg-down.csv','go.csv','pathway.csv'))file.copy(file.path(pp,file),file.path(directory,file),overwrite=TRUE)
 for(file in c('gsea.csv','gsea-ranks.csv','pathway_gsea.csv','pathway_gsea-ranks.csv'))file.copy(file.path(directory,'profiles','strict',file),file.path(directory,file),overwrite=TRUE)
 chosen<-profile_counts[[active]]
 write_json(chosen,file.path(directory,'active-profile.json'),auto_unbox=TRUE,pretty=TRUE)
 r<-read.csv(file.path(directory,'differential.csv'),check.names=FALSE)

 summary[[cmp$name]]<-list(name=cmp$name,control=cmp$control,treatment=cmp$treatment,tested=nrow(r),up=sum(r$direction=='up'),down=sum(r$direction=='down'),significance=chosen$column,alpha=chosen$alpha,lfc=chosen$lfc,shrinkage=cfg$shrinkage,missing_pvalue=sum(!is.finite(r$pvalue)),independent_filtered=sum(is.finite(r$pvalue)&!is.finite(r$padj)),rank_ties=sum(duplicated(r$stat[is.finite(r$stat)])),active_profile=active,profiles=profile_counts)
}
# Preserve first-contrast paths for existing projects and scripts.
first<-file.path(root,'contrasts',cfg$comparisons[[1]]$name)
for(f in list.files(first,full.names=TRUE)[!dir.exists(list.files(first,full.names=TRUE))])file.copy(f,file.path(root,basename(f)),overwrite=TRUE)
run_mining(vs,norm,m,mm,annotation,root,cfg)
write_json(list(design=paste(deparse(design),collapse=''),genes_input=before,genes_retained=nrow(dds),samples=nrow(m),pca_genes=length(pca_genes),reference=cfg$reference,contrasts=summary,annotation_warning=annotation$warning,gsea_engine=paste0('fgsea ',packageVersion('fgsea'),' multilevel; Wald stat; BH')),file.path(root,'analysis-summary.json'),auto_unbox=TRUE,pretty=TRUE)
writeLines(capture.output(sessionInfo()),file.path(root,'R-session.txt'));saveRDS(dds,file.path(root,'dds.rds'))
