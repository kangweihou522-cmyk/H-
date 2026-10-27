# Annotation, over-representation and fgsea outputs shared by every contrast.
load_annotation <- function(genes, root, cfg) {
 ids <- sub('\\..*$', '', genes)
 custom <- file.path(root, 'annotation.tsv')
 warning <- NULL
 if (file.exists(custom)) {
  a <- read.delim(custom, colClasses='character', na.strings='', check.names=FALSE)
  names(a)[names(a)=='gene'] <- 'id'
 } else if (requireNamespace('org.Mm.eg.db', quietly=TRUE)) {
  db <- org.Mm.eg.db::org.Mm.eg.db
  k <- intersect(ids, AnnotationDbi::keys(db, keytype='ENSEMBL'))
  if (!length(k)) stop('No Ensembl mouse IDs matched org.Mm.eg.db. Check reference / identifiers.')
  a <- AnnotationDbi::select(db, keys=k, keytype='ENSEMBL', columns=c('SYMBOL','ENTREZID','GOALL','ONTOLOGYALL'))
  a <- a[,c('ENSEMBL','SYMBOL','ENTREZID','GOALL','ONTOLOGYALL')]
  names(a) <- c('id','symbol','entrez','go','ontology')
 } else {
  a <- data.frame(id=ids, symbol=NA_character_, entrez=NA_character_, go=NA_character_, ontology=NA_character_)
  warning <- 'GO annotation unavailable: install org.Mm.eg.db or supply a real offline annotation TSV.'
 }
 a$id <- sub('\\..*$', '', a$id)
 a <- a[a$id %in% ids, , drop=FALSE]
 ann <- unique(a[!is.na(a$go) & nzchar(a$go), c('id','go','ontology')])
 if (!nrow(ann)) warning <- 'GO not run: zero genes with GO annotation.'
 basic <- a[!duplicated(a$id), , drop=FALSE]
 mapped <- data.frame(gene=genes, id=ids, symbol=basic$symbol[match(ids,basic$id)], entrez=basic$entrez[match(ids,basic$id)])
 write.csv(mapped, file.path(root,'gene-annotation.csv'), row.names=FALSE)
 go <- lapply(split(ann$id, ann$go), unique)
 descriptions <- vapply(names(go), function(term) {
  x <- GO.db::GOTERM[[term]]
  if (is.null(x)) term else AnnotationDbi::Term(x)
 }, character(1))
 ontology <- ann$ontology[match(names(go), ann$go)]; names(ontology) <- names(go)
 gmt <- list(); gd <- character()
 if (file.exists(file.path(root,'pathways.gmt'))) {
  lines <- strsplit(readLines(file.path(root,'pathways.gmt')), '\t')
  gmt <- lapply(lines, function(x) unique(x[-c(1,2)]))
  names(gmt) <- vapply(lines, `[`, character(1), 1)
  gd <- vapply(lines, `[`, character(1), 2); names(gd) <- names(gmt)
 }
 list(mapped=mapped, go=go, go_desc=descriptions, ontology=ontology, gmt=gmt, gmt_desc=gd, warning=warning)
}
empty_ora <- function() data.frame(term=character(),description=character(),ontology=character(),direction=character(),overlap=integer(),size=integer(),selected=integer(),background=integer(),gene_ratio=numeric(),rich_factor=numeric(),pvalue=numeric(),padj=numeric(),genes=character())
ora <- function(sets, desc, universe, selected, ontology, cfg) {
 # Each ontology has its own annotation coverage and hypergeometric universe.
 if(length(ontology) && length(unique(ontology))>1){
  return(do.call(rbind,lapply(sort(unique(ontology)),function(o){
   keys<-names(ontology)[ontology==o]
   ora(sets[keys],desc,universe,selected,ontology[keys],cfg)
  })))
 }
 universe <- intersect(unique(universe), unique(unlist(sets)))
 selected <- intersect(unique(selected), universe)
 eligible <- lapply(sets, intersect, y=universe)
 eligible <- eligible[lengths(eligible)>=cfg$min_set & lengths(eligible)<=cfg$max_set]
 if (!length(eligible) || !length(universe)) return(empty_ora())
 rows <- lapply(names(eligible), function(term) {
  genes <- intersect(eligible[[term]],selected); k <- length(genes); size <- length(eligible[[term]])
  data.frame(term=term,description=unname(desc[term]),ontology=if(length(ontology)) unname(ontology[term]) else 'GMT',direction='',overlap=k,size=size,selected=length(selected),background=length(universe),gene_ratio=if(length(selected)) k/length(selected) else 0,rich_factor=k/size,pvalue=phyper(k-1,size,length(universe)-size,length(selected),lower.tail=FALSE),padj=NA_real_,genes=paste(genes,collapse=';'))
 })
 out <- do.call(rbind,rows)
 out$padj <- ave(out$pvalue,out$ontology,FUN=function(x)p.adjust(x,'BH'))
 out
}
empty_gsea <- function() data.frame(term=character(),description=character(),ontology=character(),size=integer(),ES=numeric(),NES=numeric(),pvalue=numeric(),padj=numeric(),log2err=numeric(),leadingEdge=character())
gsea <- function(sets,desc,ranks,ontology,cfg) {
 if (!length(sets) || length(ranks)<cfg$min_set+1) return(empty_gsea())
 set.seed(42)
 out <- fgsea::fgseaMultilevel(pathways=sets,stats=ranks,minSize=cfg$min_set,maxSize=min(cfg$max_set,length(ranks)-1),eps=1e-10,nPermSimple=cfg$permutations,nproc=1)
 if (!nrow(out)) return(empty_gsea())
 data.frame(term=out$pathway,description=unname(desc[out$pathway]),ontology=if(length(ontology))unname(ontology[out$pathway]) else 'GMT',size=out$size,ES=out$ES,NES=out$NES,pvalue=out$pval,padj=out$padj,log2err=out$log2err,leadingEdge=vapply(out$leadingEdge,paste,character(1),collapse=';'))
}
run_enrichment <- function(r,annotation,directory,cfg,run_gsea=TRUE) {
 # ID collapse uses the highest mean expression, not the most significant statistic.
 r <- r[order(-r$baseMean,r$gene,na.last=TRUE),]
 ids <- sub('\\..*$', '',r$gene)
 map <- annotation$mapped[match(r$gene,annotation$mapped$gene),]
 gene_types <- list(go=ids, pathway=if(cfg$gmt_id=='ENSEMBL')ids else if(cfg$gmt_id=='SYMBOL')map$symbol else map$entrez)
 for (database in c('go','pathway')) {
  sets <- if(database=='go')annotation$go else annotation$gmt
  desc <- if(database=='go')annotation$go_desc else annotation$gmt_desc
  ont <- if(database=='go')annotation$ontology else character()
  mapped <- gene_types[[database]]
  valid <- !is.na(mapped) & nzchar(mapped) & is.finite(r$stat)
  keep <- valid;keep[valid] <- !duplicated(mapped[valid])
  ranks <- sort(setNames(r$stat[keep],mapped[keep]),decreasing=TRUE)
  if(database=='pathway' && length(sets) && !length(intersect(names(ranks),unique(unlist(sets)))))stop('GMT genes match zero tested genes; check gene ID type and annotation.')
  selection <- list(all=r$significant,up=r$significant & r$log2FoldChange>0,down=r$significant & r$log2FoldChange<0)
  all <- lapply(names(selection),function(direction){
   x<-ora(sets,desc,names(ranks),mapped[keep & selection[[direction]]],ont,cfg)
   if(nrow(x))x$direction<-direction
   x
  })
  out <- do.call(rbind,all);write.csv(out,file.path(directory,paste0(database,'.csv')),row.names=FALSE)
  if(!run_gsea)next
  gr <- gsea(sets,desc,ranks,ont,cfg)
  name <- if(database=='go')'gsea' else 'pathway_gsea'
  write.csv(gr,file.path(directory,paste0(name,'.csv')),row.names=FALSE)
  write.csv(data.frame(gene=names(ranks),stat=unname(ranks)),file.path(directory,paste0(name,'-ranks.csv')),row.names=FALSE)
 }
}
