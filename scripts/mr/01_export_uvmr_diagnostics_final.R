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
neutro_domain_entry('scripts/mr/01_export_uvmr_diagnostics_final.R')
# Derived candidate; see provenance/script_selection.csv. Invoke via run_job.py.
options(stringsAsFactors = FALSE)

root <- "."
mr_dir <- file.path(root, "data/gwas_cache")
supp_dir <- file.path(root, "work/mr_export")
workspace_file <- file.path(mr_dir, ".RData")
harmonised_file <- file.path(mr_dir, "harmonised_list.rds")
locked_uvmr_file <- file.path(mr_dir, "uvmr_results_full.csv")
pipeline_file <- file.path(mr_dir, "mr_noUKBB修正版.R")
output_file <- file.path(supp_dir, "Table_S4_UVMR_Diagnostics_Final.csv")
provenance_file <- file.path(supp_dir, "Table_S4_UVMR_Diagnostics_Provenance.txt")

stopifnot(neutro_file_exists(workspace_file), neutro_file_exists(harmonised_file),
          neutro_file_exists(locked_uvmr_file), neutro_file_exists(pipeline_file))

final_env <- new.env(parent = emptyenv())
loaded <- neutro_load(workspace_file, envir = final_env)
required <- c("harmonised_list", "heterogeneity_list", "pleiotropy_list",
              "steiger_list", "uvmr_results_list")
missing_objects <- setdiff(required, loaded)
if (length(missing_objects)) {
  stop("Final workspace is missing: ", paste(missing_objects, collapse = ", "))
}

retained_harmonised <- neutro_readRDS(harmonised_file)
exposure_order <- c("neutrophil", "lymphocyte", "platelet")
exposure_labels <- c(
  neutrophil = "Neutrophil count",
  lymphocyte = "Lymphocyte count",
  platelet = "Platelet count"
)
expected_ids <- c(
  neutrophil = "ieu-b-34",
  lymphocyte = "ieu-b-32",
  platelet = "ukb-d-30080_irnt"
)
expected_snps <- c(neutrophil = 383L, lymphocyte = 449L, platelet = 329L)

if (!identical(names(final_env$harmonised_list), exposure_order)) {
  stop("Unexpected exposure order in final workspace harmonised_list")
}
if (!identical(names(retained_harmonised), exposure_order)) {
  stop("Unexpected exposure order in retained harmonised_list.rds")
}

core_columns <- c(
  "SNP", "beta.exposure", "se.exposure", "effect_allele.exposure",
  "other_allele.exposure", "beta.outcome", "se.outcome",
  "effect_allele.outcome", "other_allele.outcome", "mr_keep"
)

checks <- character()
for (key in exposure_order) {
  dat_workspace <- final_env$harmonised_list[[key]]
  dat_retained <- retained_harmonised[[key]]
  if (nrow(dat_workspace) != expected_snps[[key]]) {
    stop(key, ": workspace SNP count differs from locked Table S4")
  }
  if (nrow(dat_retained) != expected_snps[[key]]) {
    stop(key, ": retained RDS SNP count differs from locked Table S4")
  }
  if (!identical(dat_workspace$SNP, dat_retained$SNP)) {
    stop(key, ": SNP vectors differ between workspace and retained RDS")
  }
  if (!isTRUE(all.equal(dat_workspace[, core_columns], dat_retained[, core_columns],
                        check.attributes = FALSE))) {
    stop(key, ": harmonized effect/allele columns differ")
  }
  if (!identical(unique(dat_workspace$id.exposure), unname(expected_ids[[key]]))) {
    stop(key, ": exposure GWAS ID differs from Methods")
  }
  if (!all(dat_workspace$mr_keep %in% TRUE)) {
    stop(key, ": final harmonized object contains mr_keep != TRUE")
  }
  if (!all(grepl("Major Depressive Disorder \\(PGC MDD2025 EUR\\)",
                 unique(dat_workspace$outcome)))) {
    stop(key, ": outcome label is not final PGC MDD2025 EUR")
  }
  checks <- c(checks, sprintf(
    "%s: %d SNPs; ID=%s; SNP/effect/allele columns identical to retained harmonised_list.rds",
    exposure_labels[[key]], nrow(dat_workspace), expected_ids[[key]]
  ))
}

pipeline_text <- paste(neutro_readLines(pipeline_file, warn = FALSE), collapse = "\n")
if (!grepl("pgc-mdd2025_no23andMe-noUKBB_eur", pipeline_text, fixed = TRUE)) {
  stop("Final UVMR pipeline does not identify the noUKBB outcome file")
}

# Verify the deterministic main estimates in the retained workspace against the
# locked UVMR CSV. Weighted-median/mode rows are intentionally not regenerated
# because their resampling implementations are stochastic and the locked values
# remain the submission source of truth.
locked_uvmr <- neutro_read_csv(locked_uvmr_file, check.names = FALSE)
method_map <- c(
  "Inverse variance weighted" = "Inverse variance weighted",
  "MR Egger" = "MR Egger"
)
for (key in exposure_order) {
  object_rows <- final_env$uvmr_results_list[[key]]
  for (method in names(method_map)) {
    obj <- object_rows[object_rows$method == method, , drop = FALSE]
    csv <- locked_uvmr[locked_uvmr$Exposure == key &
                         locked_uvmr$Method == method_map[[method]], , drop = FALSE]
    if (nrow(obj) != 1L || nrow(csv) != 1L) stop("Missing deterministic main-result row")
    if (obj$nsnp != csv$nSNP || round(obj$b, 4) != round(csv$Beta, 4) ||
        round(obj$se, 4) != round(csv$SE, 4) ||
        signif(obj$pval, 3) != signif(as.numeric(csv$Pval), 3)) {
      stop(key, " ", method, ": retained workspace main estimate conflicts with locked CSV")
    }
  }
}

rows <- list()
`%+%` <- function(a, b) paste0(a, b)
append_row <- function(...) {
  rows[[length(rows) + 1L]] <<- data.frame(..., check.names = FALSE)
}

for (key in exposure_order) {
  label <- exposure_labels[[key]]
  snps <- expected_snps[[key]]

  het <- final_env$heterogeneity_list[[key]]
  for (method in c("Inverse variance weighted", "MR Egger")) {
    h <- het[het$method == method, , drop = FALSE]
    if (nrow(h) != 1L) stop(key, ": missing heterogeneity row for ", method)
    append_row(
      Exposure = label,
      Diagnostic = if (method == "Inverse variance weighted")
        "Cochran’s Q (IVW)" else "Cochran’s Q (MR-Egger)",
      Statistic = h$Q,
      SE = NA_real_,
      df = h$Q_df,
      P = h$Q_pval,
      `Direction/Result` = NA_character_,
      SNPs = snps,
      Source_object = paste0("data/gwas_cache/.RData::heterogeneity_list[[\"", key, "\"]]")
    )
  }

  pleio <- final_env$pleiotropy_list[[key]]
  if (nrow(pleio) != 1L) stop(key, ": missing MR-Egger intercept row")
  append_row(
    Exposure = label,
    Diagnostic = "MR-Egger intercept",
    Statistic = pleio$egger_intercept,
    SE = pleio$se,
    df = NA_real_,
    P = pleio$pval,
    `Direction/Result` = NA_character_,
    SNPs = snps,
    Source_object = paste0("data/gwas_cache/.RData::pleiotropy_list[[\"", key, "\"]]")
  )

  steiger <- final_env$steiger_list[[key]]
  if (nrow(steiger) != 1L) stop(key, ": missing Steiger row")
  append_row(
    Exposure = label,
    Diagnostic = "Steiger directionality",
    Statistic = steiger$snp_r2.exposure,
    SE = NA_real_,
    df = NA_real_,
    P = steiger$steiger_pval,
    `Direction/Result` = if (isTRUE(steiger$correct_causal_direction))
      "Correct causal direction = TRUE; outcome R² = " %+% format(steiger$snp_r2.outcome, digits = 10)
    else "Correct causal direction = FALSE; outcome R² = " %+% format(steiger$snp_r2.outcome, digits = 10),
    SNPs = snps,
    Source_object = paste0("data/gwas_cache/.RData::steiger_list[[\"", key, "\"]]")
  )
}

diagnostics <- do.call(rbind, rows)
neutro_write_csv(diagnostics, output_file, row.names = FALSE, na = "")

provenance <- c(
  "Table S4 UVMR diagnostics provenance",
  paste0("Final workspace: ", workspace_file),
  paste0("Retained harmonized object: ", harmonised_file),
  paste0("Locked UVMR main-result CSV: ", locked_uvmr_file),
  paste0("Final pipeline: ", pipeline_file),
  "Outcome pipeline: PGC MDD2025 EUR excluding UK Biobank",
  "Diagnostics source: directly exported from saved final workspace objects; no diagnostic was recomputed.",
  "Instrument selection, clumping, harmonization, UVMR/MVMR estimation: not rerun.",
  checks,
  "Deterministic IVW and MR-Egger rows agree with the locked UVMR CSV at its retained precision.",
  paste0("Output: ", output_file)
)
neutro_writeLines(provenance, provenance_file, useBytes = TRUE)

neutro_cat(paste(provenance, collapse = "\n"), "\n")
neutro_cat("\nDiagnostics:\n")
print(diagnostics, row.names = FALSE)
