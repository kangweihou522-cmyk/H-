# Independent numerical checks on completed synthetic projects only.
suppressPackageStartupMessages(library(jsonlite))
a<-commandArgs(TRUE);root<-a[1]
m<-read.csv(file.path(root,'samples.csv'));v<-read.csv(file.path(root,'vst.csv'),row.names=1,check.names=FALSE)
x<-as.matrix(read.csv(file.path(root,'design-matrix.csv'),row.names=1,check.names=FALSE))
table<-read.csv(file.path(root,'gene-trait-adjusted.csv'))
for(i in seq_len(min(30,nrow(table)))){
 row<-table[i,];trait<-m[[row$trait]];fit<-lm(as.numeric(v[row$gene,m$sample])~x+trait-1)
 co<-summary(fit)$coefficients['trait',]
 df<-df.residual(fit);r<-unname(co['t value'])/sqrt(unname(co['t value'])^2+df)
 stopifnot(df==row$df,abs(r-row$correlation)<1e-10,abs(co['Pr(>|t|)']-row$pvalue)<1e-10)
}
stopifnot(max(abs(table$padj-p.adjust(table$pvalue,'BH')),na.rm=TRUE)<1e-12)
if(file.exists(file.path(root,'wgcna-module-enrichment.csv'))){
 out<-read.csv(file.path(root,'wgcna-module-enrichment.csv'))
 for(family in split(out,interaction(out$database,out$ontology,drop=TRUE))){
  stopifnot(max(abs(family$padj-p.adjust(family$pvalue,'BH')))<1e-12)
  expected<-with(family,phyper(overlap-1,size,background-size,selected,lower.tail=FALSE))
  stopifnot(max(abs(expected-family$pvalue))<1e-12)
 }
}
cat('PASS: adjusted Pearson / df / P agree with independent R lm; global BH; module hypergeometric and BH\n')
