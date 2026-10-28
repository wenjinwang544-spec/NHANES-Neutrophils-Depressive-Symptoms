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
neutro_domain_entry('scripts/nhanes/21_export_final.R')
# Derived candidate; see provenance/script_selection.csv. Invoke via run_job.py.
options(warn=1);out<-neutro_normalize('work/export');root<-getwd();log<-function(...)neutro_cat(...,'\n');csv<-function(x,n){p<-file.path(out,n);stopifnot(!neutro_file_exists(p));neutro_write_csv(x,p,row.names=FALSE);log('EXPORT',n,nrow(x))}
prov<-function(x,script,source,dest){data.frame(source_script_and_output=paste(script,source,paste0('export=',dest),sep=' | '),x,check.names=FALSE)}
miroot<-'work/imputation';miscript<-file.path(miroot,'run_sequence.R');mifile<-file.path(miroot,'Rubin_MI_CC_comparison.csv');mi<-neutro_read_csv(mifile)
stopifnot(all(mi$N_MI==33174),all(mi$N_CC==29360))
main<-subset(mi,outcome=='depression');main$FMI[main$model=='Model1']<-NA;main$lambda[main$model=='Model1']<-NA;main$RVI[main$model=='Model1']<-NA
csv(prov(main,miscript,mifile,'Table2_MI_CC.csv'),'Table2_MI_CC.csv')
h<-subset(mi,outcome=='depression'&model=='Model1');h$Holm_MI_5<-p.adjust(h$P_MI,method='holm');h$Holm_CC_5<-p.adjust(h$P_CC,method='holm');h<-h[,c('model','type','term','N_MI','N_CC','P_MI','Holm_MI_5','P_CC','Holm_CC_5','df_MI','dfcom')]
log('Model1 Holm: final Table2 MI pooled P and same-sample CC P; five tests each family; MI contrasts Barnard-Rubin df, dfcom=degf(design).')
csv(prov(h,miscript,mifile,'Model1_Holm5.csv'),'Model1_Holm5.csv')
sev<-subset(mi,outcome!='depression');sev$FMI[sev$model=='Model1']<-NA;csv(prov(sev,miscript,mifile,'Severity_MI.csv'),'Severity_MI.csv')
# Fixed SD effect scales from original Q/CI; no manually entered effects.
forest<-subset(mi,term=='log_Neutrophil');stopifnot(all(abs(forest$SD_OR_MI-exp(forest$beta_MI*.381754443))<1e-12));forest$fixed_SD<-.381754443
csv(prov(forest,miscript,mifile,'Figure2_plot_data.csv'),'Figure2_plot_data.csv')
specs<-list(
 c('Threshold_gradient_coefficients.csv','work/thresholds/run_B3_and_ALQ_audit.R','work/thresholds/B3_2_threshold_coefficients.csv'),
 c('Threshold_gradient_Holm8.csv','work/thresholds/run_B3_and_ALQ_audit.R','work/thresholds/B3_2_threshold_contrasts_Holm8.csv'),
 c('Threshold_case_counts.csv','work/thresholds/run_B3_and_ALQ_audit.R','work/thresholds/B3_1_case_counts_prevalence.csv'),
 c('Supplement_B2_sensitivity_coefficients.csv','work/symptoms_extended/run_B2_remaining.R','work/symptoms_extended/B2_2_joint_coefficients.csv'),
 c('Supplement_B2_sensitivity_contrasts.csv','work/symptoms_extended/run_B2_remaining.R','work/symptoms_extended/B2_2_t_contrasts.csv'),
 c('Supplement_B2_sensitivity_Holm.csv','work/symptoms_extended/run_B2_remaining.R','work/symptoms_extended/B2_2_holm_five_tests.csv'),
 c('Supplement_B2_neutrophil_only.csv','work/symptoms_extended/run_B2_remaining.R','work/symptoms_extended/B2_2a_neutrophil_only.csv'),
 c('Supplement_B2_same_sample_diabetes_hypertension.csv','work/symptoms_extended/run_B2_remaining.R','work/symptoms_extended/B2_2c_paired_same_sample.csv'),
 c('Supplement_B2_same_sample_education_PIR.csv','work/symptoms_extended/run_B2_remaining.R','work/symptoms_extended/B2_2d_paired_same_sample.csv'),
 c('Supplement_B2_same_sample_alcohol.csv','work/symptoms_extended/run_B2_remaining.R','work/symptoms_extended/B2_2e_paired_same_sample.csv'),
 c('Supplement_B2_sample_counts.csv','work/symptoms_extended/run_B2_remaining.R','work/symptoms_extended/B2_2_complete_case_counts.csv'),
 c('Supplement_antidepressant_coefficients.csv','work/analytic/run_D2_1.R','work/analytic/D2_1_coefficients.csv'),
 c('Supplement_antidepressant_contrasts.csv','work/analytic/run_D2_1.R','work/analytic/D2_1_t_contrasts.csv'),
 c('Supplement_antidepressant_counts.csv','work/analytic/run_D2_1.R','work/analytic/D2_1_antidep_counts.csv'),
 c('Supplement_SII_effects.csv','work/descriptive/run_E1.R','work/descriptive/E1_1_SII.csv'),
 c('Supplement_SII_constraints.csv','work/descriptive/run_E1.R','work/descriptive/E1_1_constraints.csv'),
 c('Supplement_SII_joint_coefficients.csv','work/descriptive/run_E1.R','work/descriptive/E1_1_joint_coefficients.csv'),
 c('Supplement_symptom_standardized_stacked.csv','work/symptoms_base/run_B1_2_to_6.R','work/symptoms_base/B1_4b_stacked_symptom_contrasts.csv'),
 c('Supplement_symptom_relative_scale.csv','work/symptoms_extended/run_B2_remaining.R','work/symptoms_extended/B2_3_quasiPoisson_symptom_ratios.csv'),
 c('Supplement_symptom_any_binary.csv','work/symptoms_extended/run_B2_remaining.R','work/symptoms_extended/B2_4_any_symptom_logistic.csv'),
 c('Supplement_2011_2018_coefficients.csv','work/symptoms_base/run_B1_2_to_6.R','work/symptoms_base/B1_6_joint_coefficients.csv'),
 c('Supplement_2011_2018_contrasts.csv','work/symptoms_base/run_B1_2_to_6.R','work/symptoms_base/B1_6_t_contrasts.csv'),
 c('Supplement_age20_counts.csv','work/clinical/run_D1.R','work/clinical/D1_4_adult_counts.csv'),
 c('Supplement_included_excluded.csv','work/clinical/run_D1.R','work/clinical/D1_5_included_excluded_comparison.csv'),
 c('Supplement_included_excluded_counts.csv','work/clinical/run_D1.R','work/clinical/D1_5_population_counts.csv'),
 c('Table1_weighted.csv','work/descriptive/run_E1.R','work/descriptive/E1_3_Table1.csv'),
 c('Table1_methods.csv','work/descriptive/run_E1.R','work/descriptive/E1_3_methods.csv'))
sources<-data.frame()
for(sp in specs){x<-neutro_read_csv(sp[3],check.names=FALSE);if(sp[1]=='Threshold_gradient_Holm8.csv'){stopifnot(all(abs(x$P_Holm_8-p.adjust(x$P_raw,'holm'))<1e-12))};csv(prov(x,sp[2],sp[3],sp[1]),sp[1]);sources<-rbind(sources,data.frame(script=sp[2],source=sp[3],destination=sp[1],rows=nrow(x)))}
