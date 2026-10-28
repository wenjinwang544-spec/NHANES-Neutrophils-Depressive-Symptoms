# Validated full-adult construction followed by subpopulation subset.
domain_svydesign<-function(ids=NULL,probs=NULL,strata=NULL,variables=NULL,fpc=NULL,data=NULL,nest=FALSE,...){
  dat<-as.data.frame(data);stacked<-anyDuplicated(dat$SEQN)>0;keyvars<-c('SEQN','SDMVPSU','SDMVSTRA','WTMEC_adj','WTMEC2YR');b<-full$variables
  if(stacked){stopifnot('outcome_type'%in%names(dat));lev<-levels(dat$outcome_type);blocks<-lapply(lev,function(l){q<-b;q$outcome_type<-factor(l,levels=lev);q});b<-do.call(rbind,blocks);row.names(b)<-NULL;ix<-match(paste(b$SEQN,b$outcome_type),paste(dat$SEQN,dat$outcome_type))}else ix<-match(b$SEQN,dat$SEQN)
  present<-!is.na(ix)
  for(v in intersect(c('SDMVPSU','SDMVSTRA','WTMEC_adj'),names(dat)))stopifnot(isTRUE(all.equal(as.numeric(b[[v]][present]),as.numeric(dat[[v]][ix[present]]),tolerance=1e-12)))
  for(v in setdiff(names(dat),keyvars))b[[v]]<-dat[[v]][ix]
  b$.audit_domain<-present
  if(stacked){des<-survey::svydesign(ids=~SDMVPSU+SEQN,strata=~SDMVSTRA,weights=~WTMEC_adj,data=b,nest=TRUE)}else{des<-full;des$variables<-b}
  des<-subset(des,.audit_domain);stopifnot(nrow(des$variables)==nrow(dat));des
 }
