# BEGIN MANAGED ENTRY POINT -- fails before study code without explicit workspace.
.neutro_source_root <- Sys.getenv("NEUTRO_SOURCE_ROOT")
if (!nzchar(.neutro_source_root)) {
 .entry <- grep("^--file=",commandArgs(FALSE),value=TRUE)
 if(length(.entry)!=1L)stop("Use Rscript with --workspace or the managed Python launcher")
 .cursor <- dirname(normalizePath(gsub("~+~"," ",sub("^--file=","",.entry),fixed=TRUE),mustWork=TRUE))
 repeat {if(file.exists(file.path(.cursor,"scripts/lib/safety.R"))){.neutro_source_root<-.cursor;break};.parent<-dirname(.cursor);if(.parent==.cursor)stop("Cannot locate managed source root");.cursor<-.parent}
}
base::source(file.path(.neutro_source_root,"scripts/lib/safety.R"),local=.GlobalEnv)
neutro_setup(.neutro_source_root)
if("--check-only" %in% commandArgs(TRUE)){cat("SAFE_ENTRY_CHECK_ONLY\n");quit(status=0L)}
# END MANAGED ENTRY POINT
neutro_domain_entry('scripts/nhanes/sex_interactions.R')
out<-file.path(getwd(),'work/sex');
suppressPackageStartupMessages({library(survey);library(mice)})
d<-neutro_readRDS('work/analytic/analytic_v3.rds');ds0<-neutro_readRDS('work/cohort/C2_full_adult_design.rds');imp<-neutro_readRDS('work/imputation/imputations_standardized.rds');s<-neutro_read_csv('work/imputation/imputation_scaling.csv');b<-ds0$variables;ix<-match(b$SEQN,d$SEQN);b$included<-!is.na(ix)
for(v in setdiff(names(d),c('SEQN','SDMVPSU','SDMVSTRA','WTMEC_adj','WTMEC2YR')))b[[v]]<-d[[v]][ix]
b$cycle_factor<-factor(b$SDDSRVYR);for(v in c('Neutrophil','Lymphocyte','Platelet'))b[[paste0('log_',v)]]<-log(b[[paste0(v,'_count')]])
extra<-c(Model1='age + bmi + smoking + cycle_factor + race_eth',Model2='age + bmi + smoking + cycle_factor + race_eth + log10_cotinine + alcohol4 + education + pir');pool<-function(q,u,df){p<-pool.scalar(q,u,n=df+1,k=1);se<-sqrt(p$t);c(beta=p$qbar,SE=se,df=p$df,P=2*pt(-abs(p$qbar/se),p$df),lower=p$qbar-qt(.975,p$df)*se,upper=p$qbar+qt(.975,p$df)*se)}
rows<-registry<-estimates<-data.frame();allfits<-list()
for(model in names(extra)){
 ff<-as.formula(paste('depression ~ (log_Neutrophil + log_Lymphocyte + log_Platelet)*sex +',extra[model]));fr<-as.formula(paste('depression ~ log_Neutrophil + log_Lymphocyte + log_Platelet + sex +',extra[model]));fn<-as.formula(paste('depression ~ log_Neutrophil*sex + log_Lymphocyte + log_Platelet +',extra[model]));fitlist<-list(full=list(),reduced=list(),Nonly=list())
 for(i in if(model=='Model1')1L else 1:20){filled<-complete(imp,i);for(v in intersect(s$variable,names(filled))){j<-match(v,s$variable);filled[[v]]<-filled[[v]]*s$scale[j]+s$center[j]};if(is.numeric(filled$education))filled$education<-factor(c('Less than high school','High school/GED','Some college','College or above')[as.integer(filled$education)],levels=c('Less than high school','High school/GED','Some college','College or above'))
  bb<-b;for(v in c('education','pir','alcohol4','diabetes','hypertension','log10_cotinine'))bb[[v]]<-filled[[v]][ix];des<-ds0;des$variables<-bb;des<-subset(des,included)
  for(type in names(fitlist)){f<-switch(type,full=ff,reduced=fr,Nonly=fn);ww<-character();m<-withCallingHandlers(svyglm(f,design=des,family=quasibinomial()),warning=function(w){ww<<-c(ww,conditionMessage(w));invokeRestart('muffleWarning')});stopifnot(m$converged);fitlist[[type]][[i]]<-m;registry<-rbind(registry,data.frame(model=model,type=type,imputation=i,N=nobs(m),design_df=degf(des),residual_df=m$df.residual,formula=paste(deparse(f),collapse=' '),warnings=paste(ww,collapse=' | ')));neutro_cat('FIT',model,i,type,'\n')}
 }
 if(model=='Model1')for(type in names(fitlist))fitlist[[type]]<-rep(fitlist[[type]],20)
 for(type in c('full','Nonly')){fits<-fitlist[[type]];term<-grep('log_Neutrophil:sex|sex:log_Neutrophil',names(coef(fits[[1]])),value=TRUE);q<-sapply(fits,function(m)coef(m)[term]);u<-sapply(fits,function(m)vcov(m)[term,term]);df<-min(sapply(fits,function(m)m$df.residual));pp<-pool(q,u,df)
  estimates<-rbind(estimates,data.frame(model=model,type=type,hypothesis='beta_neutrophil_x_sex = 0',N=33174,m=20,dfcom=df,t(pp),OR_ratio=exp(pp['beta']),OR_ratio_lower=exp(pp['lower']),OR_ratio_upper=exp(pp['upper'])))
 }
 j<-as.data.frame(mice::D1(as.mira(fitlist$full),as.mira(fitlist$reduced),dfcom=109)$result);rows<-rbind(rows,data.frame(model=model,test='All three cell x sex coefficients jointly zero',N=33174,m=20,dfcom=109,F=j$F.value,df1=j$df1,df2=j$df2,P=j[['P(>F)']],RIV=j$RIV))
 allfits[[model]]<-fitlist
}
neutro_write_csv(estimates,file.path(out,'interaction_hypotheses_comparison.csv'),row.names=FALSE);neutro_write_csv(rows,file.path(out,'global_D1_tests.csv'),row.names=FALSE);neutro_write_csv(registry,file.path(out,'model_registry.csv'),row.names=FALSE)
neutro_capture_output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
neutro_saveRDS(allfits$Model1,file.path(out,'Model1_validation_models.rds'))
