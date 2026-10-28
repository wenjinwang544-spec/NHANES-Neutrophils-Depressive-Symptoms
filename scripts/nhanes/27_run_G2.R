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
neutro_domain_entry('scripts/nhanes/27_run_G2.R')
# Derived candidate; see provenance/script_selection.csv. Invoke via run_job.py.
options(warn=1,digits=17)
suppressPackageStartupMessages(library(survey))
out<-'work/discrimination'
wcsv<-function(x,n){p<-file.path(out,n);stopifnot(!neutro_file_exists(p));neutro_write_csv(x,p,row.names=FALSE)}
neutro_cat('START',as.character(Sys.time()),'\n');print(sessionInfo())
seed<-20261013L;set.seed(seed);neutro_cat('BOOTSTRAP seed=',seed,' replicates=500; ordinary with-replacement PSU resampling within each stratum on complete C2 base; same draw for all paired models; refit each replicate; percentile CI.\n')
d<-neutro_readRDS('work/analytic/analytic_v3.rds');des<-neutro_readRDS('work/cohort/C2_full_adult_design.rds')
b<-des$variables;ix<-match(b$SEQN,d$SEQN);b$included<-!is.na(ix)
for(v in setdiff(names(d),c('SEQN','SDMVPSU','SDMVSTRA','WTMEC_adj')))b[[v]]<-d[[v]][ix]
b$cycle_factor<-factor(b$SDDSRVYR)
vcc<-c('log10_cotinine','alcohol4','education','pir','diabetes','hypertension','any_nsaid','any_corticosteroid','any_statin')
b$cc<-b$included&complete.cases(b[,vcc]);des$variables<-b
full<-subset(des,included);cc<-subset(des,cc);stopifnot(nrow(full$variables)==33174,nrow(cc$variables)==29360)
prev<-as.numeric(coef(svymean(~depression,full)));wcsv(data.frame(N=33174,weighted_prevalence=prev),'G2_prevalence.csv')
covs<-c(Model1='age + sex + bmi + smoking + cycle_factor + race_eth',Model2='age + sex + bmi + smoking + cycle_factor + race_eth + log10_cotinine + alcohol4 + education + pir',Model3='age + sex + bmi + smoking + cycle_factor + race_eth + log10_cotinine + alcohol4 + education + pir + diabetes + hypertension + any_nsaid + any_corticosteroid + any_statin')
auc<-function(y,s,w){o<-order(s);y<-y[o];s<-s[o];w<-w[o];g<-cumsum(c(TRUE,diff(s)!=0));p<-as.numeric(rowsum(w*y,g));n<-as.numeric(rowsum(w*(1-y),g));sum(p*(cumsum(n)-n+0.5*n))/(sum(p)*sum(n))}
neutro_cat('Weighted AUC = sum(case weight * control weight * [I(score_case>score_control)+0.5 I(tie)]) / (total case weight * total control weight); evaluated in fitting sample.\n')
psukey<-interaction(b$SDMVSTRA,b$SDMVPSU,drop=TRUE);pid<-as.integer(psukey);strata<-split(unique(data.frame(pid=pid,stratum=b$SDMVSTRA))$pid,unique(data.frame(pid=pid,stratum=b$SDMVSTRA))$stratum)
B<-500;mult<-matrix(0,length(levels(psukey)),B)
for(r in 1:B)for(ids in strata)mult[ids,r]<-tabulate(match(sample(ids,length(ids),replace=TRUE),ids),nbins=length(ids))
registry<-list();results<-list();boot<-list();fits<-list();k<-0
for(mod in names(covs)){
 dd<-if(mod=='Model1')full else cc;dv<-dd$variables;take<-match(dv$SEQN,b$SEQN);ww<-weights(dd);y<-dv$depression
 fs<-list(covariates=as.formula(paste('depression ~',covs[mod])),joint=as.formula(paste('depression ~ log_Neutrophil + log_Lymphocyte + log_Platelet +',covs[mod])))
 xs<-list();a0<-numeric(2)
 for(j in 1:2){id<-paste(mod,names(fs)[j],sep='_');neutro_cat('MODEL',id,'N',nrow(dv),'\n');print(fs[[j]])
 m<-svyglm(fs[[j]],design=dd,family=quasibinomial());fits[[id]]<-m;a0[j]<-auc(y,as.numeric(predict(m,type='response')),ww);xs[[j]]<-model.matrix(m)
 k<-k+1;registry[[k]]<-data.frame(model=id,replicate=0,N=nrow(dv),formula=paste(deparse(fs[[j]]),collapse=' '),converged=m$converged,warning='')
 }
 br<-matrix(NA_real_,B,3)
 for(r in 1:B){wr<-ww*mult[pid[take],r];keep<-wr>0
 for(j in 1:2){warnings<-character();m<-withCallingHandlers(glm.fit(xs[[j]][keep,,drop=FALSE],y[keep],weights=wr[keep]/mean(wr[keep]),family=quasibinomial()),warning=function(e){warnings<<-c(warnings,conditionMessage(e));invokeRestart('muffleWarning')})
 br[r,j]<-auc(y[keep],m$fitted.values,wr[keep]);k<-k+1;registry[[k]]<-data.frame(model=paste(mod,names(fs)[j],sep='_'),replicate=r,N=sum(keep),formula=paste(deparse(fs[[j]]),collapse=' '),converged=m$converged,warning=paste(warnings,collapse=' | '))
 }
 br[r,3]<-br[r,2]-br[r,1]
 if(r%%25==0)neutro_cat('BOOTSTRAP',mod,r,'/500',as.character(Sys.time()),'\n')
 }
 ci<-quantile(br[,3],c(.025,.975));results[[mod]]<-data.frame(model=mod,N=nrow(dv),AUC_covariates=a0[1],AUC_joint=a0[2],delta_AUC=a0[2]-a0[1],delta_lower=ci[1],delta_upper=ci[2],bootstrap_replicates=B,seed=seed)
 boot[[mod]]<-data.frame(model=mod,replicate=1:B,AUC_covariates=br[,1],AUC_joint=br[,2],delta_AUC=br[,3])
 wcsv(boot[[mod]],paste0('G2_bootstrap_',mod,'.csv'))
}
wcsv(do.call(rbind,results),'G2_1_AUC.csv');wcsv(do.call(rbind,registry),'G2_model_registry.csv');neutro_saveRDS(fits,file.path(out,'G2_models.rds'))
# E-values only from existing MI pooled CSV values, no new MI or outcome fitting.
t<-neutro_read_csv('work/export/Table2_MI_CC.csv');s<-neutro_read_csv('work/export/Severity_MI.csv')
t<-t[t$term=='log_Neutrophil'&t$type=='coefficient',];s<-s[s$term=='log_Neutrophil'&s$type=='coefficient'&s$outcome=='PHQ9_total'&s$model%in%c('Model1','Model3'),]
efun<-function(x){r<-ifelse(x<1,1/x,x);r+sqrt(r*(r-1))}
erows<-function(z,outcome){rr<-exp(z$beta_MI*0.381754443);lo<-exp(log(z$lower_MI)*0.381754443);hi<-exp(log(z$upper_MI)*0.381754443);near<-ifelse(lo<=1&hi>=1,1,ifelse(rr>=1,lo,hi));data.frame(model=z$model,outcome=outcome,SD=0.381754443,effect=rr,lower=lo,upper=hi,CI_near_null=near,E_value=efun(rr),E_value_CI=efun(near),OR_approximated_RR=outcome=='PHQ9_ge10')}
e<-rbind(erows(t,'PHQ9_ge10'),erows(s,'PHQ9_total'));wcsv(e,'G2_2_Evalues.csv')
neutro_cat('E-value formula: R + sqrt(R*(R-1)); R=max(RR,1/RR). CI containing 1 -> E-value=1. PHQ10 OR approximated RR using rare outcome prevalence=',prev,'; total score quasi-Poisson mean ratio used directly, per fixed 1 SD. Source VanderWeele and Ding 2017 doi:10.7326/M16-2607.\n')
wcsv(data.frame(seed=seed,replicates=B,SD=0.381754443,bootstrap='within-stratum PSU with replacement, full base before domain subset',CI='paired percentile'),'G2_settings.csv')
print(do.call(rbind,results));print(e);neutro_cat('COMPLETE\n')
