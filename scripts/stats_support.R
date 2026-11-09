# Shared explicit designs, estimands, durable state and adjusted associations.
atomic_rds <- function(x,path){tmp<-paste0(path,'.tmp');saveRDS(x,tmp);if(!file.rename(tmp,path))stop('Cannot save model checkpoint')}
atomic_json <- function(x,path){tmp<-paste0(path,'.tmp');jsonlite::write_json(x,tmp,auto_unbox=TRUE,pretty=TRUE,null='null');if(!file.rename(tmp,path))stop('Cannot save analysis state')}
model_formula <- function(cfg,group=TRUE,interaction=TRUE){
 variables<-vapply(cfg$covariates,`[[`,character(1),'name')
 if(group)variables<-c(variables,if(!is.null(cfg$interaction_factor)&&nzchar(cfg$interaction_factor))paste('group',cfg$interaction_factor,sep=if(interaction)' * ' else ' + ') else 'group')
 if(length(variables))reformulate(variables) else ~1
}
prepare_metadata <- function(m,cfg){
 rownames(m)<-m$sample;m$group<-relevel(factor(m$group),cfg$control)
 for(c in cfg$covariates)m[[c$name]]<-if(c$type=='numeric')as.numeric(m[[c$name]]) else factor(m[[c$name]])
 if(!is.null(cfg$interaction_factor)&&nzchar(cfg$interaction_factor))m[[cfg$interaction_factor]]<-relevel(factor(m[[cfg$interaction_factor]]),cfg$interaction_reference)
 m
}
contrast_vector <- function(m,design,cmp,cfg=NULL){
 template<-m[rep(1,if(!is.null(cmp$kind)&&cmp$kind=='interaction')4 else 2),,drop=FALSE]
 template$group<-factor(rep(c(cmp$treatment,cmp$control),nrow(template)/2),levels=levels(m$group))
 if(!is.null(cfg$interaction_factor)&&nzchar(cfg$interaction_factor)){
  f<-cfg$interaction_factor
  template[[f]]<-factor(c(rep(cmp$factor_level,2),if(nrow(template)==4)rep(cmp$reference_level,2)),levels=levels(m[[f]]))
 }
 x<-model.matrix(design,template);vector<-x[1,]-x[2,]
 if(nrow(x)==4)vector<-vector-(x[3,]-x[4,])
 vector
}
feature_traits <- function(values,traits,m,mm,cfg,key){
 rows<-list();if(!ncol(traits))return(data.frame())
 subject<-if(is.null(cfg$subject_column))'subject' else cfg$subject_column
 if(subject %in% names(m)&&length(unique(m[[subject]]))<nrow(m))return(data.frame())
 adjust<-if(is.null(cfg$trait_adjustment))'design' else cfg$trait_adjustment
 nuisance<-if(adjust=='design')mm else model.matrix(model_formula(cfg,group=FALSE),m)
 for(method in c('raw','adjusted')){
  X<-if(method=='raw')matrix(1,nrow(m),1) else nuisance
  rank<-qr(X)$rank;df<-nrow(m)-rank-1
  if(df<3)next
  projection<-diag(nrow(m))-X%*%MASS::ginv(X)
  expr<-values%*%projection
  for(trait in colnames(traits)){
   v<-traits[,trait];res<-projection%*%v
   if(sqrt(sum(res^2))<=max(1,sqrt(sum(v^2)))*1e-10)next
   den<-sqrt(rowSums(expr^2))*sqrt(sum(res^2));r<-as.vector(expr%*%res)/den;r[den<1e-10]<-NA_real_;r<-pmax(-1,pmin(1,r))
   p<-2*pt(-abs(r)*sqrt(df/pmax(1-r^2,1e-30)),df)
   z<-atanh(pmax(-1+1e-14,pmin(1-1e-14,r)));err<-qnorm(.975)/sqrt(df-1)
   out<-data.frame(feature=rownames(values),trait=trait,correlation=r,pvalue=p,method=method,n=nrow(m),df=df,ci_low=tanh(z-err),ci_high=tanh(z+err),adjustment=if(method=='raw')'Intercept only' else paste(colnames(nuisance),collapse=';'))
   names(out)[1]<-key;rows[[length(rows)+1]]<-out
  }
 }
 if(!length(rows))return(data.frame())
 out<-do.call(rbind,rows);out$padj<-ave(out$pvalue,out$method,FUN=function(x)p.adjust(x,'BH'));out
}
safe_enrichment <- function(...,directory){
 tryCatch({run_enrichment(...,directory=directory);atomic_json(list(state='completed'),file.path(directory,'enrichment-status.json'))},error=function(e){
  for(f in c('go.csv','pathway.csv'))if(!file.exists(file.path(directory,f)))write.csv(empty_ora(),file.path(directory,f),row.names=FALSE)
  atomic_json(list(state='failed',reason=conditionMessage(e)),file.path(directory,'enrichment-status.json'));warning(conditionMessage(e))
 })
}
