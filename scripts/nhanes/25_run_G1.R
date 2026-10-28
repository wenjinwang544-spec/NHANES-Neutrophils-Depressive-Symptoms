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
neutro_domain_entry('scripts/nhanes/25_run_G1.R')
# Derived candidate; see provenance/script_selection.csv. Invoke via run_job.py.
options(warn=1,digits=17)
suppressPackageStartupMessages(library(survey))
out<-neutro_normalize('work/covariates')
csv_new<-function(x,n){p<-file.path(out,n);stopifnot(!neutro_file_exists(p));neutro_write_csv(x,p,row.names=FALSE,na='NA')}
save_new<-function(x,n){p<-file.path(out,n);stopifnot(!neutro_file_exists(p));neutro_saveRDS(x,p)}
neutro_cat('START_G1',as.character(Sys.time()),'\n');print(sessionInfo())
neutro_cat('Preflight read-only R check: survey available; N=33174; complete cases=29360; cotinine nonmissing=32822; cotinine<0.015=8794; C2 base N=42143, df=109. No preflight models fitted.\n')
neutro_cat('Full C2 design then subset only. Interaction coefficient t df=degf(design); 2-df interaction joint Wald F denominator df=degf(design). Coefficient CI/P: model residual df, as D2. No additional analysis variants.\n')
inputs<-c('work/analytic/analytic_v3.rds','work/cohort/C2_full_adult_design.rds')
hash_before<-neutro_md5(inputs)
d<-as.data.frame(neutro_readRDS(inputs[1]));design<-neutro_readRDS(inputs[2]);stopifnot(nrow(d)==33174,!anyDuplicated(d$SEQN))
SDfixed<-0.381754443
for(v in c('Neutrophil','Lymphocyte','Platelet'))d[[paste0('log_',v)]]<-log(d[[paste0(v,'_count')]])
d$Neutrophil_10e9_L<-d$Neutrophil_count/1000
d$Lymphocyte_10e9_L<-d$Lymphocyte_count/1000
d$Platelet_10e9_L<-d$Platelet_count/1000
b<-design$variables;ix<-match(b$SEQN,d$SEQN);b$included<-!is.na(ix)
for(v in setdiff(names(d),c('SEQN','SDMVPSU','SDMVSTRA','WTMEC_adj')))b[[v]]<-d[[v]][ix]
b$cycle_factor<-factor(b$SDDSRVYR)
m2vars<-c('log10_cotinine','alcohol4','education','pir')
m3vars<-c(m2vars,'diabetes','hypertension','any_nsaid','any_corticosteroid','any_statin')
b$m2_complete<-b$included&complete.cases(b[,m2vars])
b$m3_complete<-b$included&complete.cases(b[,m3vars])
design$variables<-b
full<-subset(design,included);cc2<-subset(design,m2_complete);cc3<-subset(design,m3_complete)
cot<-subset(full,!is.na(LBXCOT))
stopifnot(nrow(full$variables)==33174,nrow(cc3$variables)==29360,nrow(cot$variables)==32822)
cov1<-'age + sex + bmi + smoking + cycle_factor + race_eth'
cov_without_smoking<-'age + sex + bmi + cycle_factor + race_eth'
joint_cells<-'log_Neutrophil + log_Lymphocyte + log_Platelet'
extra<-c(Model1='',Model2='+ log10_cotinine + alcohol4 + education + pir',Model3='+ log10_cotinine + alcohol4 + education + pir + diabetes + hypertension + any_nsaid + any_corticosteroid + any_statin')
registry<-data.frame();coefficients<-data.frame();fits<-list()
fit<-function(id,f,des){
 neutro_cat('\nMODEL_START',id,'N',nrow(des$variables),'\n');print(f);print(des$call)
 ww<-character();err<-''
 m<-tryCatch(withCallingHandlers(svyglm(f,design=des,family=quasibinomial()),warning=function(w){ww<<-c(ww,conditionMessage(w));neutro_cat('WARNING',id,conditionMessage(w),'\n');invokeRestart('muffleWarning')}),error=function(e){err<<-conditionMessage(e);neutro_cat('ERROR',id,err,'\n');NULL})
 if(!is.null(m)){print(summary(m));print(confint(m));fits[[id]]<<-m
  cs<-coef(summary(m));ci<-confint(m);terms<-rownames(cs)
  coefficients<<-rbind(coefficients,data.frame(model=id,term=terms,beta=unname(cs[,1]),SE=unname(cs[,2]),P=unname(cs[,ncol(cs)]),lower_beta=ci[terms,1],upper_beta=ci[terms,2]))
 }
 registry<<-rbind(registry,data.frame(model=id,N_domain=nrow(des$variables),N_fitted=if(is.null(m))NA else nobs(m),formula=paste(deparse(f),collapse=' '),family='quasibinomial',design_df=degf(des),residual_df=if(is.null(m))NA else m$df.residual,converged=if(is.null(m))NA else m$converged,warning_count=length(ww),warnings=paste(ww,collapse=' | '),error=err))
 neutro_cat('MODEL_END',id,'\n');m
}
effect<-function(m,term,delta=1){
 if(is.null(m))return(data.frame(beta=NA_real_,SE=NA_real_,OR=NA_real_,lower=NA_real_,upper=NA_real_,P=NA_real_,CI_df=NA_real_))
 beta<-unname(coef(m)[term]);se<-unname(sqrt(vcov(m)[term,term]));ci<-confint(m)[term,];cs<-coef(summary(m))
 data.frame(beta=beta,SE=se,OR=exp(beta*delta),lower=exp(ci[1]*delta),upper=exp(ci[2]*delta),P=unname(cs[term,ncol(cs)]),CI_df=m$df.residual)
}
csv_new(data.frame(sample=c('all','Model2_complete','Model3_common_complete','cotinine_nonmissing'),N=c(nrow(full$variables),nrow(cc2$variables),nrow(cc3$variables),nrow(cot$variables)),fixed_log_N_SD=SDfixed),'G1_sample_counts.csv')
# G1-1: six nested fits and four additional unique single-block fits.
# The cotinine-only model is identical in (a)/(b), and is reused, not fitted twice.
layers<-c('','+ log10_cotinine','+ log10_cotinine + alcohol4','+ log10_cotinine + alcohol4 + education + pir','+ log10_cotinine + alcohol4 + education + pir + diabetes + hypertension',extra['Model3'])
steps<-c('Model 1','+log10_cotinine','+alcohol4','+education+pir (Model 2)','+diabetes+hypertension','+three medications (Model 3)')
decomp<-data.frame()
for(i in seq_along(layers)){
 id<-paste0('G1_1a_L',i-1);m<-fit(id,as.formula(paste('depression ~',joint_cells,'+',cov1,layers[i])),cc3)
 decomp<-rbind(decomp,cbind(data.frame(section='cumulative',step=steps[i],model=id,N=nrow(cc3$variables)),effect(m,'log_Neutrophil')))
}
blocks<-c('+ log10_cotinine','+ alcohol4','+ education + pir','+ diabetes + hypertension','+ any_nsaid + any_corticosteroid + any_statin')
for(i in seq_along(blocks)){
 if(i==1){id<-'G1_1a_L1';m<-fits[[id]];neutro_cat('MODEL_REUSED G1_1b_cotinine G1_1a_L1\n')}else{id<-paste0('G1_1b_block',i);m<-fit(id,as.formula(paste('depression ~',joint_cells,'+',cov1,blocks[i])),cc3)}
 decomp<-rbind(decomp,cbind(data.frame(section='single_block',step=paste('Model 1',blocks[i]),model=id,N=nrow(cc3$variables)),effect(m,'log_Neutrophil')))
}
beta0<-decomp$beta[1];decomp$percent_change_from_Model1<-100*(decomp$beta-beta0)/beta0
csv_new(decomp,'G1_1_covariate_decomposition.csv')
# G1-2: smoking strata in all participants and those with cotinine measured.
smoking_results<-interaction_items<-interaction_joint<-data.frame()
stratum_fit<-function(id,des,sample,stratum){
 m<-fit(id,as.formula(paste('depression ~',joint_cells,'+',cov_without_smoking)),des)
 cases<-sum(des$variables$depression==1)
 smoking_results<<-rbind(smoking_results,cbind(data.frame(sample=sample,stratum=stratum,model=id,N=nrow(des$variables),PHQ10_cases=cases,fixed_SD=SDfixed,unstable_cases_lt100=cases<100,design_df=degf(des)),effect(m,'log_Neutrophil',SDfixed)))
}
for(sample in c('all','cotinine_nonmissing')){
 des<-if(sample=='all')full else cot
 for(sm in c('Never','Former','Current'))stratum_fit(paste('G1_2a',sample,sm,sep='_'),subset(des,smoking==sm),sample,sm)
 id<-paste0('G1_2a_interaction_',sample)
 m<-fit(id,as.formula(paste('depression ~',joint_cells,'+',cov1,'+ log_Neutrophil:smoking')),des)
 if(!is.null(m)){
  terms<-grep('log_Neutrophil:smoking|smoking.*:log_Neutrophil',names(coef(m)),value=TRUE)
  for(t in terms){beta<-unname(coef(m)[t]);se<-unname(sqrt(vcov(m)[t,t]));df<-degf(des);q<-qt(.975,df)
   interaction_items<-rbind(interaction_items,data.frame(sample=sample,model=id,term=t,beta=beta,SE=se,t=beta/se,df=df,lower=beta-q*se,upper=beta+q*se,P=2*pt(-abs(beta/se),df)))
  }
  wt<-regTermTest(m,~log_Neutrophil:smoking,method='Wald',df=degf(des));print(wt)
  interaction_joint<-rbind(interaction_joint,data.frame(sample=sample,N=nrow(des$variables),model=id,Wald_F=as.numeric(wt$Ftest),df1=wt$df,df2=degf(des),P=as.numeric(wt$p)))
 }else{interaction_joint<-rbind(interaction_joint,data.frame(sample=sample,N=nrow(des$variables),model=id,Wald_F=NA,df1=2,df2=degf(des),P=NA))}
}
low<-subset(full,!is.na(LBXCOT)&LBXCOT<0.015)
stratum_fit('G1_2b_cotinine_belowLOD',low,'cotinine_belowLOD','LBXCOT<0.015')
neverlow<-subset(full,smoking=='Never'&!is.na(LBXCOT)&LBXCOT<0.015)
stratum_fit('G1_2c_never_belowLOD',neverlow,'Never_and_cotinine_belowLOD','Never & LBXCOT<0.015')
csv_new(smoking_results,'G1_2_smoking_stratified.csv')
csv_new(interaction_items,'G1_2_smoking_interaction_t.csv')
csv_new(interaction_joint,'G1_2_smoking_interaction_joint_Wald.csv')
# G1-3: linear counts in 10^9/L, separate model-specific complete-case domains.
linear_results<-data.frame()
for(id in names(extra)){
 des<-if(id=='Model1')full else if(id=='Model2')cc2 else cc3
 for(type in c('neutrophil_only','joint_linear_counts')){
  exposure<-if(type=='neutrophil_only')'Neutrophil_10e9_L'else 'Neutrophil_10e9_L + Lymphocyte_10e9_L + Platelet_10e9_L'
  key<-paste('G1_3',id,type,sep='_');m<-fit(key,as.formula(paste('depression ~',exposure,'+',cov1,extra[id])),des)
  linear_results<-rbind(linear_results,cbind(data.frame(model=id,exposure_model=type,fit_id=key,N=nrow(des$variables),PHQ10_cases=sum(des$variables$depression==1),increment_10e9_L=1),effect(m,'Neutrophil_10e9_L')))
 }
}
csv_new(linear_results,'G1_3_linear_counts.csv')
csv_new(registry,'G1_model_registry.csv');csv_new(coefficients,'G1_all_model_coefficients.csv');save_new(fits,'G1_models.rds')
stopifnot(nrow(registry)==26,all(registry$N_domain==registry$N_fitted | is.na(registry$N_fitted)))
hash_after<-neutro_md5(inputs);stopifnot(all(hash_before==hash_after))
csv_new(data.frame(path=inputs,md5_before=unname(hash_before),md5_after=unname(hash_after),unchanged=hash_before==hash_after),'G1_input_integrity.csv')
csv_new(data.frame(models_run=nrow(registry),converged=sum(registry$converged,na.rm=TRUE),warning_models=sum(registry$warning_count>0),error_models=sum(nzchar(registry$error)),status='COMPLETE_WAITING_USER'),'G1_status.csv')
neutro_cat('G1_COMPLETE_WAITING_USER\n');print(decomp);print(smoking_results);print(interaction_joint);print(linear_results)
