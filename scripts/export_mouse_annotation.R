# Export org.Mm.eg.db with propagated GO memberships for offline use.
suppressPackageStartupMessages({library(org.Mm.eg.db);library(AnnotationDbi)})
a<-commandArgs(TRUE);if(length(a)!=1)stop('Usage: Rscript export_mouse_annotation.R annotation.tsv')
db<-org.Mm.eg.db;ann<-select(db,keys=keys(db,'ENSEMBL'),keytype='ENSEMBL',columns=c('SYMBOL','ENTREZID','GOALL','ONTOLOGYALL'));names(ann)<-c('gene','symbol','entrez','go','ontology')
write.table(ann,a[1],sep='\t',quote=FALSE,row.names=FALSE,na='')
