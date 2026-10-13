suppressPackageStartupMessages({library(DESeq2);library(tximport);library(GO.db);library(jsonlite)})
a <- commandArgs(TRUE); root <- a[1]; cfg <- fromJSON(file.path(root,'config.json'))
m <- read.csv(file.path(root,'samples.csv'),check.names=FALSE); rownames(m)<-m$sample
if('batch' %in% names(m)) m$batch<-factor(m$batch)
m$group <- relevel(factor(m$group),cfg$control)
design <- if('batch' %in% names(m) && length(unique(m$batch))>1) ~batch+group else ~group
if(file.exists(file.path(root,'counts.csv'))){
 c <- read.csv(file.path(root,'counts.csv'),row.names=1,check.names=FALSE)
 dds <- DESeqDataSetFromMatrix(round(as.matrix(c[,m$sample,drop=FALSE])),m,design)
}else{
 files <- file.path(root,'quant',m$sample,'quant.sf');names(files)<-m$sample
 map <- read.delim(file.path(root,'tx2gene.tsv'),header=TRUE,colClasses='character')
 txi<-tximport(files,type='salmon',tx2gene=map[,1:2],ignoreTxVersion=TRUE)
 dds<-DESeqDataSetFromTximport(txi,m,design)
}
dds<-dds[rowSums(counts(dds)>=10)>=2,]; if(nrow(dds)<10)stop('Too few expressed genes')
dds<-DESeq(dds,quiet=TRUE)
r<-as.data.frame(results(dds,contrast=c('group',cfg$treatment,cfg$control),alpha=cfg$alpha));r$gene<-rownames(r)
ids<-sub('\\..*$','',r$gene)
if(file.exists(file.path(root,'annotation.tsv'))){
 custom<-read.delim(file.path(root,'annotation.tsv'),colClasses='character')
 ann<-data.frame(ENSEMBL=custom$gene,GO=custom$go,ONTOLOGY=custom$ontology)
 r$symbol<-custom$symbol[match(ids,custom$gene)]
 map_custom<-function(col) custom[[col]][match(ids,custom$gene)]
} else if(requireNamespace('org.Mm.eg.db',quietly=TRUE)){
 db<-org.Mm.eg.db::org.Mm.eg.db
 keys_ok<-intersect(ids,keys(db,keytype='ENSEMBL'))
 if(!length(keys_ok))stop('No Ensembl mouse gene IDs matched annotation; check gene identifiers')
 syms<-mapIds(db,keys=keys_ok,column='SYMBOL',keytype='ENSEMBL',multiVals='first');r$symbol<-unname(syms[ids])
 ann<-select(db,keys=keys_ok,columns=c('GOALL','ONTOLOGYALL'),keytype='ENSEMBL');names(ann)<-c('ENSEMBL','GO','ONTOLOGY')
 map_custom<-function(col){ mapped<-mapIds(db,keys=keys_ok,column=if(col=='symbol') 'SYMBOL' else 'ENTREZID',keytype='ENSEMBL',multiVals='first');unname(mapped[ids]) }
} else {
 r$symbol<-NA_character_;ann<-data.frame(ENSEMBL=character(),GO=character(),ONTOLOGY=character())
 map_custom<-function(col) rep(NA_character_,length(ids))
 writeLines('GO annotation unavailable: install org.Mm.eg.db or provide annotation TSV. GO not run.',file.path(root,'annotation-warning.txt'))
}
ann<-ann[!is.na(ann$GO)&nzchar(ann$GO)&ann$ENSEMBL%in%ids,]
if(!nrow(ann))writeLines('GO not run: zero annotated genes.',file.path(root,'annotation-warning.txt'))
write.csv(r,file.path(root,'differential.csv'),row.names=FALSE)
vs <- varianceStabilizingTransformation(dds,blind=FALSE)
write.csv(assay(vs),file.path(root,'vst.csv'))
write.csv(counts(dds,normalized=TRUE),file.path(root,'normalized.csv'))
# GO ORA: expressed, mapped genes form the universe; each ontology has its own BH correction.

valid<-!is.na(r$padj); sig<-ids[valid & r$padj<cfg$alpha & abs(r$log2FoldChange)>=cfg$lfc]
universe<-unique(ann$ENSEMBL);sig<-intersect(sig,universe)
sets<-split(ann$ENSEMBL,ann$GO); sets<-lapply(sets,unique)
go<-lapply(names(sets),function(term){g<-intersect(sets[[term]],universe);k<-length(intersect(g,sig)); if(length(g)<5)return(NULL);data.frame(term=term,description=Term(GOTERM[[term]]),ontology=ann$ONTOLOGY[match(term,ann$GO)],overlap=k,size=length(g),pvalue=phyper(k-1,length(g),length(universe)-length(g),length(sig),lower.tail=FALSE))})
go<-do.call(rbind,go);if(is.null(go))go<-data.frame(term=character(),description=character(),ontology=character(),overlap=integer(),size=integer(),pvalue=numeric())
go$padj<-ave(go$pvalue,go$ontology,FUN=function(x)p.adjust(x,'BH'));write.csv(go,file.path(root,'go.csv'),row.names=FALSE)
# Ranked GO GSEA, weighted running sum and gene-label permutations; deterministic seed.
ranked<-r$stat;names(ranked)<-ids;ranked<-sort(ranked[is.finite(ranked)&!duplicated(ids)],decreasing=TRUE)
es<-function(x,hit){w<-abs(x); if(!any(hit)||all(hit)||sum(w[hit])==0)return(0);walk<-cumsum(ifelse(hit,w/sum(w[hit]),-1/sum(!hit)));walk[which.max(abs(walk))]}
set.seed(42);gsea<-lapply(names(sets),function(term){hit<-names(ranked)%in%sets[[term]];n<-sum(hit);if(n<15||n>500||n>=length(ranked))return(NULL);obs<-es(ranked,hit);null<-replicate(cfg$permutations,es(ranked,sample(hit)));same<-null[sign(null)==sign(obs)];nes<-if(length(same))obs/mean(abs(same)) else NA_real_;p<-(1+sum(abs(same)>=abs(obs)))/(1+length(same));data.frame(term=term,description=Term(GOTERM[[term]]),size=n,ES=obs,NES=nes,pvalue=p)})
gsea<-do.call(rbind,gsea);if(is.null(gsea))gsea<-data.frame(term=character(),description=character(),size=integer(),ES=numeric(),NES=numeric(),pvalue=numeric());gsea$padj<-p.adjust(gsea$pvalue,'BH');write.csv(gsea,file.path(root,'gsea.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(root,'R-session.txt'));saveRDS(dds,file.path(root,'dds.rds'))
if(file.exists(file.path(root,'pathways.gmt'))){
 lines<-strsplit(readLines(file.path(root,'pathways.gmt')),'\t'); ps<-lapply(lines,function(x)unique(x[-c(1,2)]));names(ps)<-sapply(lines,`[`,1);descs<-sapply(lines,`[`,2);names(descs)<-names(ps)
 mapped<-if(cfg$gmt_id=='ENSEMBL') ids else map_custom(if(cfg$gmt_id=='SYMBOL') 'symbol' else 'entrez')
 bg<-intersect(unique(mapped[!is.na(mapped)]),unique(unlist(ps)));if(!length(bg))stop('GMT genes match zero tested genes; check ID type and annotation');selected<-intersect(mapped[valid & r$padj<cfg$alpha & abs(r$log2FoldChange)>=cfg$lfc],bg)
 pr<-lapply(names(ps),function(term){g<-intersect(ps[[term]],bg);k<-length(intersect(g,selected));if(length(g)<5)return(NULL);data.frame(term=term,description=descs[term],overlap=k,size=length(g),pvalue=phyper(k-1,length(g),length(bg)-length(g),length(selected),lower.tail=FALSE))})
 pr<-do.call(rbind,pr);if(is.null(pr))pr<-data.frame(term=character(),description=character(),overlap=integer(),size=integer(),pvalue=numeric());pr$padj<-p.adjust(pr$pvalue,'BH');write.csv(pr,file.path(root,'pathway.csv'),row.names=FALSE)
 rank2<-r$stat;names(rank2)<-mapped;keep<-is.finite(rank2)&!is.na(mapped);rank2<-rank2[keep];rank2<-rank2[order(abs(rank2),decreasing=TRUE)];rank2<-sort(rank2[!duplicated(names(rank2))],decreasing=TRUE)
 set.seed(42);pg<-lapply(names(ps),function(term){hit<-names(rank2)%in%ps[[term]];n<-sum(hit);if(n<15||n>500||n>=length(rank2))return(NULL);obs<-es(rank2,hit);null<-replicate(cfg$permutations,es(rank2,sample(hit)));same<-null[sign(null)==sign(obs)];data.frame(term=term,description=descs[term],size=n,ES=obs,NES=if(length(same))obs/mean(abs(same)) else NA_real_,pvalue=(1+sum(abs(same)>=abs(obs)))/(1+length(same)))})
 pg<-do.call(rbind,pg);if(is.null(pg))pg<-data.frame(term=character(),description=character(),size=integer(),ES=numeric(),NES=numeric(),pvalue=numeric());pg$padj<-p.adjust(pg$pvalue,'BH');write.csv(pg,file.path(root,'pathway_gsea.csv'),row.names=FALSE)
}
