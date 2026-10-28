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
neutro_domain_entry('scripts/nhanes/06_prepare_supplement_inputs.R')
d<-neutro_readRDS('work/cohort/C1_rebuilt_aligned_v3.rds')
for(v in c('Neutrophil','Lymphocyte','Platelet'))d[[paste0('log_',v)]]<-log(d[[paste0(v,'_count')]])
d$race_eth<-factor(d$race_eth,levels=c('Non-Hispanic White','Hispanic','Non-Hispanic Black','Other/Multiracial'));d$education<-factor(d$education,levels=c('Less than high school','High school/GED','Some college','College or above'))
d$depression_strict<-as.numeric(d$PHQ9_total>=15)
d$SII<-d$Neutrophil_count*d$Platelet_count/d$Lymphocyte_count;d$log_SII<-log(d$SII)
neutro_saveRDS(d,'work/cohort/supplement_input.rds')
