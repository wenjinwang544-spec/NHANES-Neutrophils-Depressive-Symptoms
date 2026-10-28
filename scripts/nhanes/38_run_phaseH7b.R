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
neutro_domain_entry('scripts/nhanes/38_run_phaseH7b.R')
# Derived candidate; see provenance/script_selection.csv. Invoke via run_job.py.
options(warn=1)
suppressPackageStartupMessages({library(survey);library(mice)})
out <- 'work/sleep'; logf <- neutro_file(file.path(out,'run_phaseH7b_20261017_final4.log'),'wt')
neutro_sink(logf,split=TRUE); on.exit({neutro_sink();close(logf)},add=TRUE)
neutro_cat('START',as.character(Sys.time()),'\n')
inputs <- c('work/analytic/analytic_v3.rds','work/cohort/C2_full_adult_design.rds','work/imputation/imputations_standardized.rds','work/sleep/nhanes_h7b_merged.rds')
hash0 <- neutro_sha256(inputs); print(hash0); print(sessionInfo())
d <- neutro_readRDS(inputs[1]); base_design <- neutro_readRDS(inputs[2]); imp <- neutro_readRDS(inputs[3]); raw <- neutro_readRDS(inputs[4]); base <- base_design$variables
ix <- match(base$SEQN,d$SEQN); rawix <- match(base$SEQN,raw$SEQN)
base$included <- !is.na(ix)
# Preserve C2 full base and join analytic covariates by SEQN before any domain subset.
joinvars <- c('depression','age','sex','bmi','smoking','race_eth','education','pir','alcohol4','diabetes','hypertension','log10_cotinine','any_nsaid','any_corticosteroid','any_statin','SDDSRVYR')
for(v in joinvars) base[[v]] <- d[[v]][ix]
newraw <- c('asthma','emphysema','bronchitis','arthritis','cancer','thyroid','liver','kidney','sleep_hours')
for(v in newraw) base[[v]] <- raw[[v]][rawix]
# MCQ160G + MCQ160K form the cross-cycle emphysema/chronic bronchitis indicator.
# MCQ160O is excluded: CDC files contain it only in 2013-14, 2015-16, 2017-18.
base$emphysema_bronchitis <- ifelse(base$emphysema==1 | base$bronchitis==1,1,
  ifelse(!is.na(base$emphysema)&!is.na(base$bronchitis),0,NA_real_))
base$cycle_factor <- factor(base$SDDSRVYR)
base$PHQ10 <- as.numeric(base$depression)
cc_m3_vars <- c('alcohol4','education','pir','diabetes','hypertension','log10_cotinine','any_nsaid','any_corticosteroid','any_statin')
cc_m3 <- base$included & complete.cases(base[,cc_m3_vars])
newvars <- c('asthma','emphysema_bronchitis','sleep_hours','kidney','arthritis','cancer','thyroid','liver')
cc_new <- complete.cases(base[,newvars])
cc <- cc_m3 & cc_new
stopifnot(sum(cc_m3)==29360L)
stopifnot(sum(cc)==as.integer(sum(cc_m3 & cc_new)))
N <- sum(cc); cases <- sum(base$PHQ10[cc]==1,na.rm=TRUE)
# Specific missingness reports for full analytic sample and M3 complete-case sample.
miss <- do.call(rbind,lapply(newvars,function(v)data.frame(variable=v,N_analytic=nrow(d),missing_N=sum(is.na(raw[[v]])),missing_pct=100*mean(is.na(raw[[v]])),M3_CC_N=sum(cc_m3),M3_CC_missing_N=sum(cc_m3 & is.na(base[[v]])),M3_CC_missing_pct=100*mean(is.na(base[[v]][cc_m3])),included_in_Model4=TRUE)))
neutro_write_csv(miss,file.path(out,'variable_missingness_final4.csv'),row.names=FALSE)
neutro_write_csv(data.frame(M3_complete_case_N=sum(cc_m3),M3_complete_case_cases=sum(base$PHQ10[cc_m3]==1,na.rm=TRUE),Model4_shared_complete_case_N=N,Model4_cases=cases,excluded_for_new_covariates=sum(cc_m3)-N),file.path(out,'sample_counts_final4.csv'),row.names=FALSE)
# Approved 20-imputation object: restore original scales; observations remain unchanged.
sc <- neutro_read_csv('work/imputation/imputation_scaling.csv',check.names=FALSE)
records <- data.frame(); registry <- data.frame()
formulas <- list(Model3=PHQ10 ~ log_Neutrophil+log_Lymphocyte+log_Platelet+age+sex+bmi+smoking+cycle_factor+race_eth+log10_cotinine+alcohol4+education+pir+diabetes+hypertension+any_nsaid+any_corticosteroid+any_statin,
 Model4=PHQ10 ~ log_Neutrophil+log_Lymphocyte+log_Platelet+age+sex+bmi+smoking+cycle_factor+race_eth+log10_cotinine+alcohol4+education+pir+diabetes+hypertension+any_nsaid+any_corticosteroid+any_statin+asthma+emphysema_bronchitis+sleep_hours+kidney+arthritis+cancer+thyroid+liver)
# M3 and M4 are fitted on the same 29,024 complete cases. 20 imputations are used as confirmed; on this complete-case subset the imputed covariates are observed values.
for(i in 1:20){
  filled <- complete(imp,i)
  for(j in seq_len(nrow(sc))) if(sc$variable[j] %in% names(filled) && is.numeric(filled[[sc$variable[j]]])) filled[[sc$variable[j]]] <- filled[[sc$variable[j]]]*sc$scale[j]+sc$center[j]
  bb <- base
  for(v in joinvars) if(v %in% names(filled)) bb[[v]][!is.na(ix)] <- filled[[v]][match(bb$SEQN[!is.na(ix)],d$SEQN)]
  # Factor variables are sourced from the original complete-case analytic data; no replacement levels are introduced.
  bb$cycle_factor <- factor(bb$SDDSRVYR)
  # Restore unchanged values on complete cases from original analytic dataset to avoid numerical round-trip drift.
  for(v in c('depression','age','sex','bmi','smoking','race_eth','education','pir','alcohol4','diabetes','hypertension','log10_cotinine','any_nsaid','any_corticosteroid','any_statin')) bb[[v]][!is.na(ix)] <- d[[v]][match(bb$SEQN[!is.na(ix)],d$SEQN)]
  bb$PHQ10 <- as.numeric(bb$depression); bb$cycle_factor <- factor(bb$SDDSRVYR)
  des <- base_design; des$variables <- bb; des$variables$cc_h7b <- cc
  subdes <- subset(des,cc_h7b)
  stopifnot(nrow(subdes$variables)==N)
  for(id in names(formulas)){
    warnings <- character(); f <- formulas[[id]]
    neutro_cat('\nMODEL_START',id,'imputation',i,'N',N,'cases',cases,'formula',paste(deparse(f),collapse=' '),'\n')
    fit <- withCallingHandlers(svyglm(f,design=subdes,family=quasibinomial()),warning=function(w){warnings<<-c(warnings,conditionMessage(w));neutro_cat('WARNING',conditionMessage(w),'\n');invokeRestart('muffleWarning')})
    if(!isTRUE(fit$converged)||length(warnings)) stop(paste('Model warning/nonconvergence:',id,i,paste(warnings,collapse=';')))
    sm <- summary(fit); co <- coef(fit); vv <- vcov(fit); df <- degf(subdes); crit <- qt(.975,df)
    terms <- if(id=='Model3') c('log_Neutrophil','log_Lymphocyte','log_Platelet') else names(co)
    for(term in terms){b<-unname(co[term]);se<-sqrt(unname(vv[term,term]));p<-2*pt(-abs(b/se),df);records<-rbind(records,data.frame(model=id,imputation=i,term=term,beta=b,SE=se,df=df,P=p,OR=exp(b),lower=exp(b-crit*se),upper=exp(b+crit*se),N=N,cases=cases))}
    registry<-rbind(registry,data.frame(model=id,imputation=i,N=N,cases=cases,df=df,converged=fit$converged,warning_count=length(warnings),warnings=paste(warnings,collapse=' | '),formula=paste(deparse(f),collapse=' ')))
    neutro_cat('MODEL_END',id,'imputation',i,'\n')
  }
}
neutro_write_csv(records,file.path(out,'per_imputation_coefficients_final4.csv'),row.names=FALSE)
neutro_write_csv(registry,file.path(out,'model_registry_final4.csv'),row.names=FALSE)
# Rubin pooling. Since all model covariates are complete in this shared sample, imputations reproduce the same fit; pool nevertheless follows the approved m=20 workflow.
pooled <- data.frame(); keys <- unique(records[,c('model','term')])
for(kk in seq_len(nrow(keys))){
  key <- keys[kk,]; z<-records[records$model==key$model & records$term==key$term,]; q<-z$beta; u<-z$SE^2; dfcom<-unique(z$df); stopifnot(length(dfcom)==1)
  ps<-mice::pool.scalar(q,u,n=dfcom+1,k=1)
  se<-sqrt(ps$t); df<-ps$df; crit<-qt(.975,df); p<-2*pt(-abs(ps$qbar/se),df); term<-key$term
  pooled<-rbind(pooled,data.frame(model=key$model,term=term,beta=ps$qbar,SE=se,df=df,P=p,OR=exp(ps$qbar),lower=exp(ps$qbar-crit*se),upper=exp(ps$qbar+crit*se),per_SD_OR=if(term=='log_Neutrophil')exp(ps$qbar*.381754443)else NA_real_,per_SD_lower=if(term=='log_Neutrophil')exp((ps$qbar-crit*se)*.381754443)else NA_real_,per_SD_upper=if(term=='log_Neutrophil')exp((ps$qbar+crit*se)*.381754443)else NA_real_,N=N,cases=cases,m=20,between_imputation_variance=ps$b,within_variance=ps$ubar,FMI=ps$fmi))
}
neutro_write_csv(pooled,file.path(out,'Model3_Model4_pooled_results_final4.csv'),row.names=FALSE)
# Restricted per-SD comparison table.
neutro_write_csv(pooled[pooled$term=='log_Neutrophil',],file.path(out,'neutrophil_per_SD_comparison_final4.csv'),row.names=FALSE)
stopifnot(all(hash0==neutro_sha256(inputs)))
neutro_write_csv(data.frame(path=inputs,sha256_before=unname(hash0),sha256_after=unname(neutro_sha256(inputs))),file.path(out,'source_integrity_final4.csv'),row.names=FALSE)
neutro_write_csv(data.frame(added_coefficients=length(newvars),added_terms=paste(newvars,collapse='; '),physical_activity_included=FALSE,COPD_MCq160O_all_cycle=FALSE),file.path(out,'model4_definition_final4.csv'),row.names=FALSE)
neutro_cat('\nCOMPLETED',as.character(Sys.time()),'N',N,'cases',cases,'added_coefficients',length(newvars),'\n');print(pooled)
