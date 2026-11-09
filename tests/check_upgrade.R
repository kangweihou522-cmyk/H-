# Independent official-function and lm oracles on synthetic completed jobs.
suppressPackageStartupMessages({library(DESeq2);library(jsonlite)})
a<-commandArgs(TRUE);root<-a[1];kind<-a[2];cfg<-fromJSON(file.path(root,'config.json'),simplifyVector=FALSE);dds<-readRDS(file.path(root,'dds.rds'));m<-as.data.frame(colData(dds))
files<-lapply(cfg$comparisons,function(cmp)read.csv(file.path(root,'contrasts',cmp$name,'differential.csv')))
if(kind=='interaction'){
 official<-results(dds,name=tail(resultsNames(dds),1),alpha=.05)
 cmp<-which(vapply(cfg$comparisons,function(c)identical(c$kind,'interaction'),logical(1)))
 table<-files[[cmp]];match<-match(table$gene,rownames(official))
 stopifnot(max(abs(table$stat-official$stat[match]),na.rm=TRUE)<1e-8,max(abs(table$log2FoldChange-official$log2FoldChange[match]),na.rm=TRUE)<1e-8)
 stopifnot(sum(table$padj[table$gene%in%rownames(dds)[1:100]]<.05,na.rm=TRUE)>80)
 reduced<-~batch+group+sex
}else{
 for(i in seq_along(cfg$comparisons)){
  cmp<-cfg$comparisons[[i]];official<-results(dds,contrast=c('group',cmp$treatment,cmp$control),alpha=.05);table<-files[[i]];idx<-match(table$gene,rownames(official))
  stopifnot(max(abs(table$stat-official$stat[idx]),na.rm=TRUE)<1e-8,max(abs(table$log2FoldChange_raw-official$log2FoldChange[idx]),na.rm=TRUE)<1e-8)
 }
 reduced<-~batch
 for(pair in list(c('pathway-activity.csv','pathway-trait.csv','term'),c('wgcna-eigengenes.csv','wgcna-module-trait.csv','module'))){
  if(!file.exists(file.path(root,pair[2])))next
  expr<-as.matrix(read.csv(file.path(root,pair[1]),row.names=1,check.names=FALSE));if(pair[3]=='module')expr<-t(expr)
  result<-read.csv(file.path(root,pair[2]));result<-result[result$trait=='disease_score'&result$method=='adjusted',]
  mm<-as.matrix(read.csv(file.path(root,'design-matrix.csv'),row.names=1,check.names=FALSE));trait<-as.numeric(as.character(m$disease_score))
  for(i in seq_len(nrow(result))){r<-result[i,];y<-as.numeric(expr[r[[pair[3]]],m$sample]);model<-lm(y~mm+trait-1);co<-summary(model)$coefficients['trait',];df<-df.residual(model);rho<-co['t value']/sqrt(co['t value']^2+df)
   stopifnot(abs(rho-r$correlation)<1e-8,abs(co['Pr(>|t|)']-r$pvalue)<1e-8,df==r$df)
  }
 }
}
lrt<-nbinomLRT(dds,reduced=reduced,quiet=TRUE);official<-results(lrt,alpha=.05);table<-read.csv(file.path(root,'omnibus.csv'));idx<-match(table$gene,rownames(official));stopifnot(max(abs(table$stat-official$stat[idx]),na.rm=TRUE)<1e-8,max(abs(table$padj-official$padj[idx]),na.rm=TRUE)<1e-8)
p<-unlist(lapply(files,`[[`,'pvalue'));q<-unlist(lapply(files,`[[`,'padj_global'));stopifnot(max(abs(q-p.adjust(p,'BH')),na.rm=TRUE)<1e-10)
cat('PASS: official DESeq2 Wald, LRT, full-family BH; independent lm adjusted association oracle; ',kind,'\n')
