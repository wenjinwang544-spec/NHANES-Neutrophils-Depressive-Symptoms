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
neutro_domain_entry('scripts/nhanes/02_build_full_adult_design.R')
suppressPackageStartupMessages(library(survey))
a<-neutro_readRDS('work/cohort/C1_all_DEMO_reconstructed_v3.rds');b<-a[a$age>=18,];b$smoking<-factor(c('Current','Former','Never')[b$smoking],levels=c('Never','Current','Former'))
b$analysis_included<-with(b,!is.na(PHQ9_total)&!is.na(Neutrophil_count)&!is.na(Lymphocyte_count)&!is.na(Platelet_count)&Neutrophil_count>0&Lymphocyte_count>0&Platelet_count>0&!is.na(bmi)&!is.na(smoking))
for(v in c('Neutrophil','Lymphocyte','Platelet'))b[[paste0('log_',v)]]<-log(b[[paste0(v,'_count')]])
full<-svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~WTMEC_adj,data=b,nest=TRUE)
stopifnot(nrow(b)==42143,sum(b$analysis_included)==33174,degf(full)==109)
neutro_saveRDS(b,'work/cohort/C2_adult_base.rds');neutro_saveRDS(full,'work/cohort/C2_full_adult_design.rds');neutro_write_csv(data.frame(C1_passed=TRUE),'work/cohort/phase_gate_status_v3.csv',row.names=FALSE)
