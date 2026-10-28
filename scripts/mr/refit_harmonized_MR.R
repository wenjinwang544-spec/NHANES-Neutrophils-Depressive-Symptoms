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
neutro_domain_entry('scripts/mr/refit_harmonized_MR.R')
set.seed(20261027L) # Validation seed; historical MR stochastic seed was not retained.
suppressPackageStartupMessages(library(TwoSampleMR))
out<-'work/mr_refit';neutro_dir_create(out,recursive=TRUE,showWarnings=FALSE)
h<-neutro_readRDS('data/gwas_cache/harmonised_list.rds');rows<-list();diag<-list()
for(n in names(h)){d<-subset(h[[n]],mr_keep);stopifnot(nrow(d)>0);r<-mr(d,method_list=c('mr_ivw','mr_egger_regression','mr_weighted_median','mr_weighted_mode'));r$exposure_key<-n;rows[[n]]<-r;diag[[n]]<-list(heterogeneity=mr_heterogeneity(d),pleiotropy=mr_pleiotropy_test(d))}
neutro_write_csv(do.call(rbind,rows),file.path(out,'uvmr_refitted.csv'),row.names=FALSE);neutro_saveRDS(diag,file.path(out,'diagnostics.rds'))
m<-neutro_readRDS('data/gwas_cache/mvmr_harmonised_final.rds');neutro_write_csv(mv_multiple(m)$result,file.path(out,'mvmr_refitted.csv'),row.names=FALSE)
neutro_capture_output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
