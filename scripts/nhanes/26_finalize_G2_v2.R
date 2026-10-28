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
neutro_domain_entry('scripts/nhanes/26_finalize_G2_v2.R')
# Derived candidate; see provenance/script_selection.csv. Invoke via run_job.py.
out<-'work/discrimination'
r<-neutro_read_csv(file.path(out,'G2_model_registry.csv'),na.strings=character(),colClasses='character');stopifnot(all(r$converged=="TRUE"),all(r$warning==''));z<-data.frame(original_models=6,bootstrap_refits=3000,converged=sum(r$converged=="TRUE"),warnings=sum(nzchar(r$warning)),bootstrap_replicates_per_model=500);neutro_write_csv(z,file.path(out,'G2_status.csv'),row.names=FALSE);print(z)
neutro_writeLines(c('Weighted AUC: weighted case-control concordance; tied predictions receive 0.5. Evaluated within fitting sample (apparent AUC).','Bootstrap: 500 ordinary within-stratum PSU samples with replacement, same number of PSUs as original stratum; sampled from complete C2 base before analysis subset. Original MEC weights multiplied by PSU draw multiplicity. Paired models refitted each replicate; delta AUC percentile CI. Seed 20261013.','Bootstrap glm.fit estimates checked against original svyglm coefficients; validation fits logged in validation.log.','E-value: R + sqrt(R*(R-1)), R=max(effect,1/effect). CI crossing 1 gives confidence-limit E-value 1. Endpoint nearest 1 separately recorded in final CSV.','PHQ9>=10: OR approximated RR; measured weighted prevalence 0.076424070451672, below 15% rare-outcome criterion. Approximation, not exact conversion.','PHQ9_total: quasi-Poisson conditional mean ratio used directly; E-value applied on ratio scale. All effects per fixed log neutrophil SD 0.381754443.','VanderWeele TJ, Ding P. Ann Intern Med. 2017;167:268-274. doi:10.7326/M16-2607. https://pubmed.ncbi.nlm.nih.gov/28693043/','Implementation reference https://louisahsmith.github.io/evalue/reference/evalues.OR.html'),file.path(out,'G2_methods.txt'))
