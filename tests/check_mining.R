# Run from project root: docker run --rm -v "$PWD:/project:ro" mouse-rnaseq:2 Rscript /project/tests/check_mining.R
suppressPackageStartupMessages(library(jsonlite))
source('/project/scripts/mining.R')
check_skip <- function(n,subjects=NULL,reason){
 root<-tempfile();dir.create(root)
 m<-data.frame(sample=paste0('S',seq_len(n)),group=factor(rep(c('C','T'),length.out=n)))
 if(!is.null(subjects))m$subject<-subjects
 cfg<-list(pathway_activity=FALSE,wgcna=TRUE,traits=list())
 run_mining(matrix(1,1,n),matrix(1,1,n),m,model.matrix(~group,m),list(gmt=list()),root,cfg)
 s<-read_json(file.path(root,'mining-status.json'))
 stopifnot(s$wgcna$state=='skipped',grepl(reason,s$wgcna$reason))
 unlink(root,recursive=TRUE)
}
check_skip(6,reason='at least 15')
check_skip(36,rep(paste0('person',1:18),each=2),reason='Repeated subjects')
cat('PASS: six independent samples skipped; 36 repeated measurements from 18 subjects also skipped\n')
