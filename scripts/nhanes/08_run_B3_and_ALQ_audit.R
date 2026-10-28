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
neutro_domain_entry('scripts/nhanes/08_run_B3_and_ALQ_audit.R')
# Derived candidate; see provenance/script_selection.csv. Invoke via run_job.py.
options(warn = 1)
suppressPackageStartupMessages({library(survey); library(haven)})
out <- neutro_normalize("work/thresholds")
input <- neutro_normalize("work/cohort/supplement_input.rds")
csv <- function(x,n) {p<-file.path(out,n);stopifnot(!neutro_file_exists(p));neutro_write_csv(x,p,row.names=FALSE,na="")}
neutro_cat("START ",format(Sys.time(),tz="Asia/Shanghai"),"\nINPUT ",input,"\nMD5 ",unname(neutro_md5(input)),"\n")
print(sessionInfo());print(options()[c("contrasts","survey.lonely.psu")])
neutro_cat("EXPLORATORY_SEVERITY_GRADIENT thresholds=5,10,15,20; models=7; contrasts=t,df=degf(design); 8 threshold contrasts one Holm family.\n")
d<-neutro_readRDS(input);d$cycle_factor<-factor(d$SDDSRVYR,levels=4:10)
covars<-"age + sex + bmi + smoking + cycle_factor + race_eth"
cells<-c("log_Neutrophil","log_Lymphocyte","log_Platelet")
stopifnot(nrow(d)==33174,!anyNA(d$PHQ9_total),all(d$PHQ9_total>=0 & d$PHQ9_total<=27),all(d$PHQ9_total==d$PHQ9_SO+d$PHQ9_CA))
for(k in c(5,10,15,20))d[[paste0("case",k)]]<-as.numeric(d$PHQ9_total>=k)
stopifnot(all(d$case10==d$depression),all(d$case15==d$depression_strict))
ds<-svydesign(id=~SDMVPSU,strata=~SDMVSTRA,weights=~WTMEC_adj,nest=TRUE,data=d)
cases<-do.call(rbind,lapply(c(5,10,15,20),function(k){x<-d[[paste0("case",k)]]==1;p<-svymean(as.formula(paste0("~case",k)),ds);data.frame(threshold=k,N=nrow(d),cases=sum(x),weighted_prevalence=as.numeric(coef(p)),prevalence_SE=as.numeric(SE(p)),PSUs_with_cases=nrow(unique(d[x,c("SDMVSTRA","SDMVPSU")])),strata_with_cases=length(unique(d$SDMVSTRA[x])),total_nested_PSUs=nrow(unique(d[,c("SDMVSTRA","SDMVPSU")])),total_strata=length(unique(d$SDMVSTRA)))}))
csv(cases,"B3_1_case_counts_prevalence.csv")
registry<-data.frame();warnings_all<-data.frame(model=character(),warning=character())
fit<-function(id,response,design,family){
 f<-as.formula(paste(response,"~",paste(cells,collapse=" + "),"+",covars));warn<-character()
 neutro_cat("MODEL_START ",id," N=",nrow(design$variables)," formula=",paste(deparse(f),collapse=" ")," family=",family$family," link=",family$link,"\n",sep="");print(design$call)
 m<-withCallingHandlers(svyglm(f,design=design,family=family),warning=function(w){warn<<-c(warn,conditionMessage(w));neutro_cat("MODEL_WARNING ",id," ",conditionMessage(w),"\n");invokeRestart("muffleWarning")})
 print(summary(m));print(confint(m))
 registry<<-rbind(registry,data.frame(model=id,N=nobs(m),formula=paste(deparse(f),collapse=" "),family=family$family,link=family$link,df_residual=m$df.residual,df_design=degf(design),converged=isTRUE(m$converged),iterations=m$iter,warning_count=length(warn)))
 if(length(warn))warnings_all<<-rbind(warnings_all,data.frame(model=id,warning=warn))
 neutro_cat("MODEL_END ",id," CONVERGED=",isTRUE(m$converged)," WARNING_COUNT=",length(warn),"\n",sep="");m
}
extract<-function(m,id){cs<-coef(summary(m))[cells,,drop=FALSE];ci<-confint(m)[cells,,drop=FALSE];data.frame(model=id,N=nobs(m),term=cells,beta=unname(coef(m)[cells]),SE=unname(sqrt(diag(vcov(m)))[cells]),ratio=unname(exp(coef(m)[cells])),CI_lower=unname(exp(ci[,1])),CI_upper=unname(exp(ci[,2])),P_raw=unname(cs[,ncol(cs)]),df_residual=m$df.residual)}
contrasts<-function(m,id,design){one<-function(expr,label){o<-svycontrast(m,expr);b<-as.numeric(coef(o));se<-as.numeric(SE(o));t<-b/se;df<-degf(design);data.frame(model=id,N=nobs(m),comparison=label,log_difference=b,SE=se,t=t,df=df,P_raw=2*pt(-abs(t),df))};r<-rbind(one(quote(log_Neutrophil-log_Lymphocyte),"Neutrophil-Lymphocyte"),one(quote(log_Neutrophil-log_Platelet),"Neutrophil-Platelet"));print(r);r}
threshold_coef<-threshold_ct<-data.frame()
for(k in c(5,10,15,20)){id<-paste0("B3_2_PHQ",k);m<-fit(id,paste0("case",k),ds,quasibinomial());cf<-extract(m,id);ct<-contrasts(m,id,ds);cf$threshold<-k;ct$threshold<-k;threshold_coef<-rbind(threshold_coef,cf);threshold_ct<-rbind(threshold_ct,ct)}
threshold_ct$P_Holm_8<-p.adjust(threshold_ct$P_raw,"holm")
csv(threshold_coef,"B3_2_threshold_coefficients.csv");csv(threshold_ct,"B3_2_threshold_contrasts_Holm8.csv")
m<-fit("B3_3_PHQ_total_quasiPoisson","PHQ9_total",ds,quasipoisson(link="log"))
csv(extract(m,"B3_3_PHQ_total_quasiPoisson"),"B3_3_total_score_rate_ratios.csv");csv(contrasts(m,"B3_3_PHQ_total_quasiPoisson",ds),"B3_3_total_score_t_contrasts.csv")
sens<-list(A=subset(ds,Lymphocyte_count<=5000),B=subset(ds,Neutrophil_count<=9600));scf<-sct<-data.frame()
for(id in names(sens)){name<-paste0("B3_4_",id);m<-fit(name,"case15",sens[[id]],quasibinomial());scf<-rbind(scf,extract(m,name));sct<-rbind(sct,contrasts(m,name,sens[[id]]))}
csv(scf,"B3_4_PHQ15_sensitivity_coefficients.csv");csv(sct,"B3_4_PHQ15_sensitivity_contrasts.csv")
csv(registry,"B3_model_registry.csv");csv(warnings_all,"B3_model_warnings.csv")
