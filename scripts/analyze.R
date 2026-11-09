suppressPackageStartupMessages({library(DESeq2);library(tximport);library(GO.db);library(jsonlite);library(fgsea)})
a <- commandArgs(TRUE);root<-a[1];cfg<-fromJSON(file.path(root,'config.json'),simplifyVector=FALSE)
if(is.null(cfg$deg_profile))cfg$deg_profile<-'strict'
script <- sub('^--file=','',grep('^--file=',commandArgs(),value=TRUE))
# Rscript encodes spaces as ~+~ in --file. The app supplies the exact mounted path.
script_dir<-Sys.getenv('RNASEQ_SCRIPT_DIR',unset=dirname(gsub('~+~',' ',script,fixed=TRUE)))
source(file.path(script_dir,'stats_support.R'));source(file.path(script_dir,'enrichment.R'));source(file.path(script_dir,'mining.R'))
m<-read.csv(file.path(root,'samples.csv'),check.names=FALSE,colClasses='character');rownames(m)<-m$sample
m<-prepare_metadata(m,cfg)
design<-model_formula(cfg)
mm<-model.matrix(design,m);if(qr(mm)$rank<ncol(mm)||nrow(mm)<=ncol(mm))stop('Design is confounded or lacks residual degrees of freedom.')
if(file.exists(file.path(root,'counts.csv'))){
 counts<-read.csv(file.path(root,'counts.csv'),row.names=1,check.names=FALSE,colClasses=c('character',rep('numeric',nrow(m))))
 dds<-DESeqDataSetFromMatrix(as.matrix(counts[,m$sample,drop=FALSE]),m,design)
}else{
 files<-file.path(root,'quant',m$sample,'quant.sf');names(files)<-m$sample
 map<-read.delim(file.path(root,'tx2gene.tsv'),header=TRUE,colClasses='character')
 map[[1]]<-sub('\\..*$','',sub('\\|.*$','',map[[1]]))
 qnames<-read.delim(files[[1]],colClasses=c(Name='character'))$Name;normalized_ids<-sub('\\..*$','',sub('\\|.*$','',qnames));atomic_json(list(transcripts=length(qnames),matched=sum(normalized_ids%in%map[[1]]),fraction=mean(normalized_ids%in%map[[1]])),file.path(root,'tx2gene-coverage.json'))
 txi<-tximport(files,type='salmon',tx2gene=map[,1:2],ignoreTxVersion=TRUE,ignoreAfterBar=TRUE)
 write.csv(txi$abundance,file.path(root,'tpm.csv'))
 dds<-DESeqDataSetFromTximport(txi,m,design)
}
raw_all<-counts(dds);before<-nrow(dds);dds<-dds[rowSums(counts(dds)>=cfg$min_count)>=cfg$min_samples,];if(nrow(dds)<10)stop('Too few expressed genes after filtering')
checkpoint<-file.path(root,'dds.rds');marker<-file.path(root,'fit-signature.rds')
inputs<-c(file.path(root,'samples.csv'),if(file.exists(file.path(root,'counts.csv')))file.path(root,'counts.csv') else c(file.path(root,'tx2gene.tsv'),file.path(root,'quant',m$sample,'quant.sf')))
signature<-list(inputs=unname(tools::md5sum(inputs)),model=cfg[c('control','covariates','design_type','interaction_factor','interaction_reference','min_count','min_samples')],code=unname(tools::md5sum(c(file.path(script_dir,'analyze.R'),file.path(script_dir,'stats_support.R')))),R=R.version.string,DESeq2=as.character(packageVersion('DESeq2')),tximport=as.character(packageVersion('tximport')))
reused<-file.exists(checkpoint)&&file.exists(marker)&&identical(readRDS(marker),signature)
if(reused){dds<-readRDS(checkpoint);message('Reusing verified DESeq2 fit checkpoint')}else{dds<-DESeq(dds,quiet=TRUE);atomic_rds(dds,checkpoint);atomic_rds(signature,marker)}
atomic_json(list(state='completed',stage='DESeq2 fit',reused=reused,method='Input content hashes, model, script and R/package versions checked; checkpoint before optional modules'),file.path(root,'statistics-checkpoint.json'))
write.csv(mm,file.path(root,'design-matrix.csv'));write.csv(model.matrix(model_formula(cfg,group=FALSE),m),file.path(root,'trait-design-matrix.csv'))
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
summary<-list();single_de<-identical(cfg$analysis_stage,'differential')
for(cmp in cfg$comparisons){
 directory<-file.path(root,'contrasts',cmp$name);dir.create(directory,recursive=TRUE,showWarnings=FALSE)
 result<-if(is.null(cfg$interaction_factor)||!nzchar(cfg$interaction_factor))results(dds,contrast=c('group',cmp$treatment,cmp$control),alpha=cfg$alpha) else results(dds,contrast=as.numeric(contrast_vector(m,design,cmp,cfg)),alpha=cfg$alpha)
 if(isTRUE(cfg$effect_test)){vec<-if(is.null(cfg$interaction_factor)||!nzchar(cfg$interaction_factor))c('group',cmp$treatment,cmp$control) else as.numeric(contrast_vector(m,design,cmp,cfg));er<-as.data.frame(results(dds,contrast=vec,alpha=.05,lfcThreshold=1,altHypothesis='greaterAbs'));er$gene<-rownames(er);er$null_hypothesis<-'abs(true log2FC) <= 1';write.csv(er,file.path(directory,'effect-threshold.csv'),row.names=FALSE)}
 r<-as.data.frame(result);r$gene<-rownames(r);r$log2FoldChange_raw<-r$log2FoldChange;r$lfcSE_raw<-r$lfcSE
 r$ci_low_raw<-r$log2FoldChange-qnorm(.975)*r$lfcSE;r$ci_high_raw<-r$log2FoldChange+qnorm(.975)*r$lfcSE
 r$missing_reason<-ifelse(!is.finite(r$pvalue),'Cook outlier or unavailable test',ifelse(!is.finite(r$padj),'Independent filtering','Tested'))
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
  if(!single_de)safe_enrichment(r=rp,annotation=annotation,cfg=cfg,run_gsea=profile=='strict',directory=pp)
  profile_counts[[profile]]<-list(profile=profile,column=spec$column,alpha=spec$alpha,lfc=spec$lfc,operator='strict < / >',total=sum(rp$significant),up=sum(rp$direction=='up'),down=sum(rp$direction=='down'))
 }
 active<-cfg$deg_profile
 pp<-file.path(directory,'profiles',active)
 for(file in c('differential.csv','deg-all.csv','deg-up.csv','deg-down.csv','go.csv','pathway.csv'))if(file.exists(file.path(pp,file)))file.copy(file.path(pp,file),file.path(directory,file),overwrite=TRUE)
 for(file in c('gsea.csv','gsea-ranks.csv','pathway_gsea.csv','pathway_gsea-ranks.csv'))if(file.exists(file.path(directory,'profiles','strict',file)))file.copy(file.path(directory,'profiles','strict',file),file.path(directory,file),overwrite=TRUE)
 chosen<-profile_counts[[active]]
 write_json(chosen,file.path(directory,'active-profile.json'),auto_unbox=TRUE,pretty=TRUE)
 r<-read.csv(file.path(directory,'differential.csv'),check.names=FALSE)

 summary[[cmp$name]]<-list(name=cmp$name,control=cmp$control,treatment=cmp$treatment,tested=nrow(r),up=sum(r$direction=='up'),down=sum(r$direction=='down'),significance=chosen$column,alpha=chosen$alpha,lfc=chosen$lfc,shrinkage=cfg$shrinkage,missing_pvalue=sum(!is.finite(r$pvalue)),independent_filtered=sum(is.finite(r$pvalue)&!is.finite(r$padj)),rank_ties=sum(duplicated(r$stat[is.finite(r$stat)])),active_profile=active,profiles=profile_counts,kind=if(is.null(cmp$kind))'pairwise' else cmp$kind,factor_level=cmp$factor_level,reference_level=cmp$reference_level,fdr_family='Genes within this contrast',lfc_source=cfg$shrinkage)
 atomic_json(list(design=paste(deparse(design),collapse=''),samples=nrow(m),genes_input=before,genes_retained=nrow(dds),pca_genes=length(pca_genes),contrasts=summary),file.path(root,'analysis-summary.json'))
}
if(isTRUE(cfg$omnibus)){
 reduced<-if(!is.null(cfg$interaction_factor)&&nzchar(cfg$interaction_factor))model_formula(cfg,interaction=FALSE) else model_formula(cfg,group=FALSE)
 lrt<-tryCatch({object<-nbinomLRT(dds,reduced=reduced,quiet=TRUE);x<-as.data.frame(results(object,alpha=.05));x$gene<-rownames(x);x$log2FoldChange<-NULL;x$lfcSE<-NULL;x$test<-if(!is.null(cfg$interaction_factor)&&nzchar(cfg$interaction_factor))'Interaction LRT' else 'Group omnibus LRT';write.csv(x,file.path(root,'omnibus.csv'),row.names=FALSE);list(state='completed',full=paste(deparse(design),collapse=''),reduced=paste(deparse(reduced),collapse=''),direction='Omnibus has no single LFC or direction')},error=function(e)list(state='failed',reason=conditionMessage(e)))
 atomic_json(lrt,file.path(root,'omnibus-status.json'))
}
if(isTRUE(cfg$global_fdr)){
 files<-file.path(root,'contrasts',vapply(cfg$comparisons,`[[`,character(1),'name'),'differential.csv');rows<-lapply(files,read.csv);pvals<-unlist(lapply(rows,`[[`,'pvalue'));q<-p.adjust(pvals,'BH');offset<-0
 for(i in seq_along(files)){x<-rows[[i]];x$padj_global<-q[offset+seq_len(nrow(x))];offset<-offset+nrow(x);write.csv(x,files[i],row.names=FALSE)
  for(profile in c('nominal','strict')){f<-file.path(dirname(files[i]),'profiles',profile,'differential.csv');y<-read.csv(f);y$padj_global<-x$padj_global[match(y$gene,x$gene)];write.csv(y,f,row.names=FALSE)
   for(direction in c('all','up','down')){g<-file.path(dirname(f),paste0('deg-',direction,'.csv'));z<-read.csv(g);z$padj_global<-x$padj_global[match(z$gene,x$gene)];write.csv(z,g,row.names=FALSE)}
  }
 }
}
# Preserve first-contrast paths for existing projects and scripts.
first<-file.path(root,'contrasts',cfg$comparisons[[1]]$name)
for(f in list.files(first,full.names=TRUE)[!dir.exists(list.files(first,full.names=TRUE))])file.copy(f,file.path(root,basename(f)),overwrite=TRUE)
if(single_de){
 write_json(list(pathway_activity=list(state='skipped',reason='Standalone differential expression'),regulatory_activity=list(state='skipped',reason='Standalone differential expression'),wgcna=list(state='skipped',reason='Standalone differential expression')),file.path(root,'mining-status.json'),auto_unbox=TRUE,pretty=TRUE)
}else run_mining(vs,norm,m,mm,annotation,root,cfg)
write_json(list(design=paste(deparse(design),collapse=''),genes_input=before,genes_retained=nrow(dds),samples=nrow(m),pca_genes=length(pca_genes),reference=cfg$reference,contrasts=summary,annotation_warning=annotation$warning,gsea_engine=if(single_de)'Skipped in standalone differential' else paste0('fgsea ',packageVersion('fgsea'),' multilevel; Wald stat; BH')),file.path(root,'analysis-summary.json'),auto_unbox=TRUE,pretty=TRUE)
writeLines(capture.output(sessionInfo()),file.path(root,'R-session.txt'));saveRDS(dds,file.path(root,'dds.rds'))
