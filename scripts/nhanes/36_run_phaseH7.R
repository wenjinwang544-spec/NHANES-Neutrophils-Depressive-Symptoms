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
neutro_domain_entry('scripts/nhanes/36_run_phaseH7.R')
# Derived candidate; see provenance/script_selection.csv. Invoke via run_job.py.
options(warn = 1)
suppressPackageStartupMessages({library(survey); library(mice)})
out <- neutro_normalize('work/race')
logfile <- file.path(out, 'run_log.txt')
logcon <- neutro_file(logfile, open='wt')
neutro_sink(logcon, split = TRUE); neutro_sink(logcon, type = 'message')
on.exit({neutro_sink(type='message'); neutro_sink(); close(logcon)}, add=TRUE)
neutro_cat('H7 start:', format(Sys.time(), '%Y-%m-%d %H:%M:%S %Z'), '\n')
neutro_cat('Script: run_phaseH7.R\n')
inputs <- c('work/imputation/imputations_standardized.rds',
            'work/imputation/imputation_scaling.csv',
            'work/cohort/C2_full_adult_design.rds',
            'work/analytic/analytic_v3.rds')
stopifnot(all(neutro_file_exists(inputs)))
hash_before <- neutro_md5(inputs)
imp <- neutro_readRDS(inputs[1]); scaling <- neutro_read_csv(inputs[2], stringsAsFactors=FALSE)
ds0 <- neutro_readRDS(inputs[3]); dat <- as.data.frame(neutro_readRDS(inputs[4]))
stopifnot(imp$m == 20L, nrow(imp$data) == 33174L, nrow(dat) == 33174L)
fixed_sd <- 0.381754443
neutro_cat('MI source:', inputs[1], 'm=', imp$m, '; fixed SD=', fixed_sd, '\n')
neutro_cat('Survey base rows:', nrow(ds0$variables), '; C2 full design retained; analysis domain via subset.\n')
# Reconstitute the exact H7 analysis domain on the original C2 design.
base <- ds0$variables
ix <- match(base$SEQN, dat$SEQN)
base$included <- !is.na(ix)
for (v in c('race_eth','sex','smoking','alcohol4','log10_cotinine','education','pir','diabetes','hypertension',
            'any_nsaid','any_corticosteroid','any_statin','depression','log_Neutrophil','log_Lymphocyte','log_Platelet')) {
  if (v %in% names(dat)) base[[v]] <- dat[[v]][ix]
}
base$cycle_factor <- factor(base$SDDSRVYR)
base$race_eth <- factor(base$race_eth, levels=c('Non-Hispanic White','Non-Hispanic Black','Hispanic','Other/Multiracial'))
ds0$variables <- base
fulldesign <- subset(ds0, included)
stopifnot(nrow(fulldesign$variables) == 33174L)
# D2's approved inverse transformation of standardized numeric variables.
inverse_data <- function(z) {
  for (j in seq_len(nrow(scaling))) {
    v <- scaling$variable[j]
    if (v %in% names(z) && is.numeric(z[[v]])) z[[v]] <- z[[v]] * scaling$scale[j] + scaling$center[j]
  }
  if ('education' %in% names(z)) z$education <- factor(as.character(z$education),
    levels=c('Less than high school','High school/GED','Some college','College or above'))
  if ('depression' %in% names(z)) z$depression <- as.numeric(z$depression > 0.5)
  z
}
models <- list(
  Model1=list(extra=''),
  Model3=list(extra=' + log10_cotinine + alcohol4 + education + pir + diabetes + hypertension + any_nsaid + any_corticosteroid + any_statin')
)
fitlog <- data.frame(); coef_rows <- data.frame(); int_rows <- data.frame()
fit_one <- function(des, f, model_name, imputation, subset_name, outcome='depression') {
  warns <- character()
  neutro_cat('MODEL_START name=', model_name, ' subset=', subset_name, ' imputation=', imputation,
      ' n=', nrow(des$variables), ' formula=', paste(deparse(f), collapse=' '), '\n', sep='')
  fit <- withCallingHandlers(svyglm(f, design=des, family=quasibinomial(link='logit')),
    warning=function(w) { warns <<- c(warns, conditionMessage(w)); neutro_cat('MODEL_WARNING ', conditionMessage(w), '\n', sep=''); invokeRestart('muffleWarning') })
  neutro_cat('MODEL_END converged=', isTRUE(fit$converged), ' warnings=', length(warns),
      ' df_design=', degf(des), '\n', sep='')
  fitlog <<- rbind(fitlog, data.frame(model=model_name, outcome=outcome, subset=subset_name, imputation=imputation,
    N=nobs(fit), formula=paste(deparse(f),collapse=' '), converged=isTRUE(fit$converged),
    design_df=degf(des), warning_count=length(warns), warnings=paste(warns,collapse=' | ')))
  fit
}
# Fit 20 race-specific models per model; retain coefficients/variances for Rubin pooling.
race_names <- levels(fulldesign$variables$race_eth)
fits_strata <- list()
for (mi in seq_len(20)) {
  filled <- inverse_data(complete(imp, mi))
  bb <- ds0$variables
  ix0 <- match(bb$SEQN, dat$SEQN)
  bb$included <- !is.na(ix0)
  for (v in names(filled)) if (v %in% names(bb)) bb[[v]][!is.na(ix0)] <- filled[[v]][ix0[!is.na(ix0)]]
  bb$cycle_factor <- factor(bb$SDDSRVYR)
  bb$race_eth <- factor(bb$race_eth, levels=race_names)
  dsi <- ds0; dsi$variables <- bb
  domain <- subset(dsi, included)
  observed_race_counts <- table(domain$variables$race_eth)
  expected_race_counts <- table(factor(dat$race_eth, levels=race_names))
  if (!identical(as.integer(observed_race_counts), as.integer(expected_race_counts))) stop('Race group N mismatch after SEQN alignment')
  for (mn in names(models)) {
    ex <- models[[mn]]$extra
    f <- as.formula(paste('depression ~ log_Neutrophil + log_Lymphocyte + log_Platelet + age + sex + bmi + smoking + cycle_factor', ex))
    for (rr in race_names) {
      subd <- subset(domain, race_eth == rr)
      fit <- fit_one(subd, f, mn, mi, rr)
      b <- unname(coef(fit)['log_Neutrophil']); v <- unname(vcov(fit)['log_Neutrophil','log_Neutrophil'])
      fits_strata[[paste(mn,rr,mi,sep='|')]] <- c(beta=b,var=v,df=degf(subd),N=nobs(fit))
      coef_rows <- rbind(coef_rows, data.frame(model=mn,race_eth=rr,imputation=mi,N=nobs(fit),
        cases=sum(subd$variables$depression==1,na.rm=TRUE),beta=b,variance=v,dfcom=degf(subd)))
    }
  }
  # Full-sample race interaction and nested reduced model for MI Wald D1 tests.
  for (mn in names(models)) {
    ex <- models[[mn]]$extra
    rhs <- paste('log_Neutrophil * race_eth + log_Lymphocyte + log_Platelet + age + sex + bmi + smoking + cycle_factor', ex)
    f1 <- as.formula(paste('depression ~', rhs))
    f0 <- as.formula(paste('depression ~ log_Neutrophil + race_eth + log_Lymphocyte + log_Platelet + age + sex + bmi + smoking + cycle_factor', ex))
    fits_strata[[paste('INTERACTION_FULL',mn,mi,sep='|')]] <- fit_one(domain,f1,paste0(mn,'_race_interaction_full'),mi,'all')
    fits_strata[[paste('INTERACTION_NULL',mn,mi,sep='|')]] <- fit_one(domain,f0,paste0(mn,'_race_interaction_null'),mi,'all')
  }
}
# Pool scalar stratum coefficients with Barnard-Rubin df using each stratum's design df.
stratum_results <- data.frame()
for (mn in names(models)) for (rr in race_names) {
  z <- coef_rows[coef_rows$model==mn & coef_rows$race_eth==rr,]
  p <- pool.scalar(Q=z$beta,U=z$variance,n=unique(z$dfcom)[1]+1,k=1)
  se <- sqrt(p$t); crit <- qt(.975,p$df)
  stratum_results <- rbind(stratum_results,data.frame(model=mn,race_eth=rr,N=unique(z$N)[1],cases=unique(z$cases)[1],
    beta=p$qbar,SE=se,df=p$df,OR_per_fixed_SD=exp(p$qbar*fixed_sd),
    CI_lower=exp((p$qbar-crit*se)*fixed_sd),CI_upper=exp((p$qbar+crit*se)*fixed_sd),P=2*pt(-abs(p$qbar/se),p$df),
    fixed_SD=fixed_sd,MI=20))
}
neutro_write_csv(stratum_results,file.path(out,'H7_1_race_stratified_MI.csv'),row.names=FALSE)
# Rubin D1 joint Wald test for the three neutrophil-by-race interaction coefficients.
interaction_results <- data.frame()
for (mn in names(models)) {
  ff <- lapply(seq_len(20), function(i) fits_strata[[paste('INTERACTION_FULL',mn,i,sep='|')]])
  fr <- lapply(seq_len(20), function(i) fits_strata[[paste('INTERACTION_NULL',mn,i,sep='|')]])
  d1 <- mice::D1(mice::as.mira(ff), mice::as.mira(fr), dfcom=degf(fulldesign))
  print(d1)
  vv <- as.numeric(d1$result[1,])
  interaction_results <- rbind(interaction_results,data.frame(model=mn,test='joint log_Neutrophil x race_eth (3 df)',
    F=vv[1],df1=vv[2],df2=vv[3],P=vv[4],RIV=vv[5],MI=20,design_df=degf(fulldesign)))
}
neutro_write_csv(interaction_results,file.path(out,'H7_1_race_interaction_MI.csv'),row.names=FALSE)
neutro_write_csv(fitlog,file.path(out,'model_registry.csv'),row.names=FALSE)
# Local availability audit; do not fit a partial model and label it Model 4.
extra_audit <- data.frame(
  requested=c('Sleep duration (hours/night)','Physical activity','Self-reported asthma',
    'Self-reported COPD/emphysema/chronic bronchitis','Self-reported arthritis','Self-reported cancer',
    'Self-reported chronic kidney disease (weak/failing kidneys)','Self-reported thyroid disease','Self-reported liver disease'),
  source_fields=c('SLD010H (2005-06 to 2013-14); SLD012 (2015-16 and 2017-18; weekday/workday); SLD013 (2017-18 weekend)',
    'No PAQ activity measure found in available C1 variable-label inventory', 'MCQ010 (Ever been told you have asthma)',
    'MCQ160G (emphysema) + MCQ160K (chronic bronchitis); combined COPD flag not present',
    'MCQ160A (doctor ever said you had arthritis)', 'MCQ220 (ever told cancer/malignancy)',
    'No weak/failing-kidney self-report variable found in available C1 variable-label inventory',
    'MCQ160M (ever told thyroid problem)', 'MCQ160L (ever told any liver condition)'),
  available_cycles=c('2005-06 to 2017-18 (field names differ)', 'not found in available C1 variable-label inventory',
    rep('2005-06 to 2017-18',1),'2005-06 to 2017-18', '2005-06 to 2017-18', '2005-06 to 2017-18',
    'not found in available C1 variable-label inventory','2005-06 to 2017-18','2005-06 to 2017-18'),
  definition=c('reported sleep hours; 2017-18 separate weekday/workday and weekend fields',
    'no harmonizable measure available in local analytic files', 'ever told by a doctor/health professional had asthma',
    'ever told had emphysema and/or chronic bronchitis; no direct combined COPD item in local label inventory',
    'doctor ever said arthritis', 'ever told cancer or malignancy',
    'weak/failing kidneys definition unavailable; no matching variable found locally',
    'ever told thyroid problem', 'ever told any liver condition'),
  missingness_in_analytic_v3=rep('not computable: requested field absent from analytic_v3.rds',9),
  status=c('sleep fields present in analytic_v3 but not in existing 20 MI data; raw special values require harmonization',
    'excluded: no local measure/data', 'not in analytic_v3/C2/MI; source XPT for 2005-18 unavailable locally',
    'not in analytic_v3/C2/MI; source XPT for 2005-18 unavailable locally',
    'not in analytic_v3/C2/MI; source XPT for 2005-18 unavailable locally',
    'not in analytic_v3/C2/MI; source XPT for 2005-18 unavailable locally',
    'excluded: no matching self-report variable found locally',
    'not in analytic_v3/C2/MI; source XPT for 2005-18 unavailable locally',
    'not in analytic_v3/C2/MI; source XPT for 2005-18 unavailable locally'),
  stringsAsFactors=FALSE)
# Replace sleep audit counts with observed RDS counts, clearly limited to raw NA coding.
sleep_cycle <- ifelse(dat$SDDSRVYR %in% 4:8,dat$SLD010H,ifelse(dat$SDDSRVYR %in% 9:10,dat$SLD012,NA))
extra_audit$missingness_in_analytic_v3[1] <- sprintf('%d/%d (%.4f%%) NA only; special values not recoded',sum(is.na(sleep_cycle)),length(sleep_cycle),100*mean(is.na(sleep_cycle)))
neutro_write_csv(extra_audit,file.path(out,'H7_2_extra_covariate_availability.csv'),row.names=FALSE)
# Cycle-specific sleep missingness; available named raw fields and NA-only rate.
cyc <- sort(unique(dat$SDDSRVYR)); cycmap <- c('4'='2005-06','5'='2007-08','6'='2009-10','7'='2011-12','8'='2013-14','9'='2015-16','10'='2017-18')
sleep_rows <- do.call(rbind,lapply(cyc,function(cy){s <- if(cy<=8) dat$SLD010H[dat$SDDSRVYR==cy] else dat$SLD012[dat$SDDSRVYR==cy];data.frame(cycle=cycmap[as.character(cy)],variable=if(cy<=8)'SLD010H' else 'SLD012',N=length(s),missing_N=sum(is.na(s)),missing_pct=100*mean(is.na(s)))}))
neutro_write_csv(sleep_rows,file.path(out,'H7_2_sleep_by_cycle.csv'),row.names=FALSE)
# Model 4 was not estimated: listed covariates are unavailable in source analysis data.
neutro_write_csv(data.frame(model='Model 4',status='NOT_FIT',reason='Required expanded covariates are absent from analytic_v3.rds, C2 full design and 20 MI object; 2005-2018 source XPT files are not locally available. Fitting M3+sleep alone would not be the prespecified Model 4.',N=NA_integer_),file.path(out,'H7_2_model4_status.csv'),row.names=FALSE)
neutro_cat('\nH7 extra variable audit:\n'); print(extra_audit); print(sleep_rows)
neutro_cat('Input MD5 before/after:\n'); print(data.frame(path=inputs,before=unname(hash_before),after=unname(neutro_md5(inputs))))
stopifnot(identical(unname(hash_before),unname(neutro_md5(inputs))))
neutro_cat('H7 end:',format(Sys.time(),'%Y-%m-%d %H:%M:%S %Z'),'\n')
