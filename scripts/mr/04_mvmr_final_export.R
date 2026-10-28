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
neutro_domain_entry('scripts/mr/04_mvmr_final_export.R')
# Derived candidate; see provenance/script_selection.csv. Invoke via run_job.py.
library(TwoSampleMR)
library(dplyr)
library(MVMR)

input_dir <- "data/gwas_cache"
out_dir <- "work/mr_export"
workspace_file <- file.path(input_dir, ".RData")
result_file <- file.path(input_dir, "mvmr_results.csv")

exposure_ids <- c("ieu-b-34", "ieu-b-32", "ukb-d-30080_irnt")
exposure_labels <- c(
  "ieu-b-34" = "Neutrophil",
  "ieu-b-32" = "Lymphocyte",
  "ukb-d-30080_irnt" = "Platelet"
)

# Load the exact harmonized MVMR object retained from the final analysis run.
final_env <- new.env(parent = emptyenv())
neutro_load(workspace_file, envir = final_env)
if (!exists("mvmr_harmonised", envir = final_env, inherits = FALSE)) {
  stop("The final workspace does not contain mvmr_harmonised")
}
mvmr_harmonised <- final_env$mvmr_harmonised

# mvmr_results.csv is the sole source for reported coefficient estimates.
saved_results <- neutro_read_csv(result_file, stringsAsFactors = FALSE)

# Recompute only as an integrity check that the saved result and retained
# harmonized object belong to the same final MVMR analysis.
check_res <- mv_multiple(mvmr_harmonised)$result
check_columns <- c("id.exposure", "b", "se", "pval", "nsnp")
saved_check <- saved_results[order(saved_results$id.exposure), check_columns]
recomputed_check <- check_res[order(check_res$id.exposure), check_columns]
rownames(saved_check) <- NULL
rownames(recomputed_check) <- NULL
if (!isTRUE(all.equal(saved_check, recomputed_check, tolerance = 1e-12))) {
  stop("mvmr_results.csv does not match mv_multiple(mvmr_harmonised)")
}

mvmr_format <- format_mvmr(
  BXGs = mvmr_harmonised$exposure_beta,
  BYG = mvmr_harmonised$outcome_beta,
  seBXGs = mvmr_harmonised$exposure_se,
  seBYG = mvmr_harmonised$outcome_se,
  RSID = rownames(mvmr_harmonised$exposure_beta)
)
mvmr_strength <- strength_mvmr(r_input = mvmr_format, gencov = 0)

strength_df <- data.frame(
  id.exposure = colnames(mvmr_harmonised$exposure_beta),
  Conditional_F = as.numeric(mvmr_strength[1, ]),
  stringsAsFactors = FALSE
)

total_snps <- nrow(mvmr_harmonised$exposure_beta)

final_mvmr <- saved_results %>%
  transmute(
    id.exposure,
    Exposure = unname(exposure_labels[id.exposure]),
    Beta = b,
    SE = se,
    OR = exp(b),
    CI_lower = exp(b - 1.96 * se),
    CI_upper = exp(b + 1.96 * se),
    P_value = pval,
    Exposure_specific_instruments = nsnp,
    SNPs_used_in_joint_model = total_snps
  ) %>%
  left_join(strength_df, by = "id.exposure") %>%
  mutate(id.exposure = factor(id.exposure, levels = exposure_ids)) %>%
  arrange(id.exposure) %>%
  mutate(id.exposure = as.character(id.exposure))

stopifnot(
  nrow(final_mvmr) == 3,
  identical(final_mvmr$id.exposure, exposure_ids),
  all(!is.na(final_mvmr$Conditional_F)),
  all(final_mvmr$SNPs_used_in_joint_model == total_snps)
)

neutro_write_csv(
  final_mvmr,
  file.path(out_dir, "Supplementary_Table_MVMR_Results.csv"),
  row.names = FALSE
)
neutro_saveRDS(
  mvmr_harmonised,
  file.path(out_dir, "mvmr_harmonised_final.rds")
)

neutro_cat("mvmr_results.csv matches the retained final mvmr_harmonised object.\n")
neutro_cat("\nFinal MVMR table:\n")
print(final_mvmr)
neutro_cat("\nTotal SNPs used in the joint MVMR model:", total_snps, "\n")
