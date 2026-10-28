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
neutro_domain_entry('scripts/nhanes/sex_total_score.R')
out<-file.path(getwd(),'work/sex_total');
suppressPackageStartupMessages({library(survey);library(mice)});saved<-neutro_readRDS('work/sex/Model1_validation_models.rds');des<-saved$full[[1]]$survey.design;stopifnot(nrow(des$variables)==33174);extra<-'age + bmi + smoking + cycle_factor + race_eth'
ff<-as.formula(paste('PHQ9_total ~ (log_Neutrophil + log_Lymphocyte + log_Platelet)*sex +',extra));fr<-as.formula(paste('PHQ9_total ~ log_Neutrophil + log_Lymphocyte + log_Platelet + sex +',extra));fn<-as.formula(paste('PHQ9_total ~ log_Neutrophil*sex + log_Lymphocyte + log_Platelet +',extra));fits<-lapply(list(full=ff,reduced=fr,Nonly=fn),function(f)svyglm(f,design=des,family=quasipoisson(link='log')));rows<-data.frame()
for(type in c('full','Nonly')){m<-fits[[type]];t<-grep('log_Neutrophil:sex',names(coef(m)),value=TRUE);q<-unname(coef(m)[t]);u<-unname(vcov(m)[t,t]);p<-pool.scalar(rep(q,20),rep(u,20),n=m$df.residual+1,k=1);se<-sqrt(p$t);rows<-rbind(rows,data.frame(model='Model1',outcome='PHQ9_total',type=type,N=nobs(m),beta=p$qbar,SE=se,dfcom=m$df.residual,df=p$df,P=2*pt(-abs(p$qbar/se),p$df),lower=p$qbar-qt(.975,p$df)*se,upper=p$qbar+qt(.975,p$df)*se,archived_design_df_P=2*pt(-abs(q/sqrt(u)),109)))}
j<-as.data.frame(D1(as.mira(rep(list(fits$full),20)),as.mira(rep(list(fits$reduced),20)),dfcom=109)$result);neutro_write_csv(rows,file.path(out,'total_score_interactions.csv'),row.names=FALSE);neutro_write_csv(j,file.path(out,'total_score_global_D1.csv'),row.names=FALSE);print(rows);print(j)
