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
neutro_domain_entry('scripts/nhanes/07_run_B2_remaining.R')
# Derived candidate; see provenance/script_selection.csv. Invoke via run_job.py.
options(warn = 1)
suppressPackageStartupMessages(library(survey))
out <- neutro_normalize("work/symptoms_extended")
input <- neutro_normalize("work/cohort/supplement_input.rds")
csv <- function(x, name) {
  p <- file.path(out, name); if (neutro_file_exists(p)) stop("Output exists: ", p)
  neutro_write_csv(x, p, row.names = FALSE, na = "")
}
neutro_cat("START ", format(Sys.time(), tz = "Asia/Shanghai"), "\nINPUT ", input, "\nMD5 ", unname(neutro_md5(input)), "\n")
print(sessionInfo()); print(options()[c("contrasts", "survey.lonely.psu", "survey.adjust.domain.lonely")])
neutro_cat("Contrasts=t; df=degf(design); coefficients=summary.svyglm residual-df t; CI=confint.svyglm.\n")
d <- neutro_readRDS(input); d$cycle_factor <- factor(d$SDDSRVYR, levels = 4:10)
d$SO_any <- as.numeric(d$PHQ9_SO >= 1); d$CA_any <- as.numeric(d$PHQ9_CA >= 1)
basevars <- c("depression", "log_Neutrophil", "log_Lymphocyte", "log_Platelet", "age", "sex", "bmi", "smoking", "cycle_factor", "race_eth", "WTMEC_adj", "SDMVPSU", "SDMVSTRA", "SEQN")
stopifnot(nrow(d) == 33174L, !anyNA(d[basevars]), sum(d$depression == 1) == 2902L, levels(d$smoking)[1] == "Never")
ds <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC_adj, nest = TRUE, data = d)
covars <- "age + sex + bmi + smoking + cycle_factor + race_eth"
cells <- c("log_Neutrophil", "log_Lymphocyte", "log_Platelet")
formula_joint <- function(extra = "", outcome = "depression") as.formula(paste(outcome, "~", paste(cells, collapse = " + "), "+", covars, extra))
registry <- data.frame(); warnings_all <- data.frame(model = character(), warning = character())
fit <- function(id, formula, design, family) {
  warnings <- character()
  neutro_cat("MODEL_START ", id, " N=", nrow(design$variables), " formula=", paste(deparse(formula), collapse = " "), " family=", family$family, " link=", family$link, "\n", sep = "")
  print(design$call)
  m <- withCallingHandlers(svyglm(formula, design = design, family = family), warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); neutro_cat("MODEL_WARNING ", id, " ", conditionMessage(w), "\n"); invokeRestart("muffleWarning") })
  print(summary(m)); print(confint(m))
  registry <<- rbind(registry, data.frame(model = id, formula = paste(deparse(formula), collapse = " "), N = nobs(m), family = family$family, link = family$link, df_residual = m$df.residual, df_design = degf(design), converged = isTRUE(m$converged), iterations = m$iter, warning_count = length(warnings)))
  if (length(warnings)) warnings_all <<- rbind(warnings_all, data.frame(model = id, warning = warnings))
  neutro_cat("MODEL_END ", id, " CONVERGED=", isTRUE(m$converged), " ITERATIONS=", m$iter, " WARNING_COUNT=", length(warnings), "\n", sep = "")
  m
}
extract <- function(m, terms, id) {
  cs <- coef(summary(m))[terms, , drop = FALSE]; ci <- confint(m)[terms, , drop = FALSE]
  data.frame(model = id, N = nobs(m), term = terms, beta = unname(coef(m)[terms]), SE = unname(sqrt(diag(vcov(m)))[terms]), OR = unname(exp(coef(m)[terms])), CI_lower = unname(exp(ci[, 1])), CI_upper = unname(exp(ci[, 2])), P_raw = unname(cs[, ncol(cs)]), df_residual = m$df.residual)
}
results <- list()
run_joint <- function(id, des, extra = "", outcome = "depression") {
  m <- fit(id, formula_joint(extra, outcome), des, quasibinomial())
  cf <- extract(m, cells, id)
  contrast <- function(expr, label) {
    obj <- svycontrast(m, expr); b <- as.numeric(coef(obj)); se <- as.numeric(SE(obj)); t <- b / se; df <- degf(des)
    data.frame(model = id, N = nobs(m), comparison = label, beta_difference = b, SE = se, t = t, df = df, P_raw = 2 * pt(-abs(t), df), CI_lower = b - qt(.975, df) * se, CI_upper = b + qt(.975, df) * se)
  }
  ct <- rbind(contrast(quote(log_Neutrophil - log_Lymphocyte), "Neutrophil-Lymphocyte"), contrast(quote(log_Neutrophil - log_Platelet), "Neutrophil-Platelet"))
  hm <- data.frame(model = id, test = c(cells, ct$comparison), P_raw = c(cf$P_raw, ct$P_raw)); hm$P_Holm <- p.adjust(hm$P_raw, "holm")
  cf$P_Holm <- hm$P_Holm[1:3]; ct$P_Holm <- hm$P_Holm[4:5]
  neutro_cat("SVYCONTRAST ", id, " t df=", degf(des), "\n"); print(ct); print(hm)
  results[[id]] <<- list(coefficients = cf, contrasts = ct, holm = hm)
}
# B2-1 M3 output is reused; no additional full-sample M3 model is fitted.
priorcf <- neutro_read_csv(file.path(out, "B2_1_M3_joint_coefficients.csv"))
priorct <- neutro_read_csv(file.path(out, "B2_1_M3_t_contrasts.csv"))
priorhm <- neutro_read_csv(file.path(out, "B2_1_M3_holm_five_tests.csv"))
priorcf$model <- "M3_reference"; priorhm$model <- "M3_reference"
priorcf$P_Holm <- priorhm$P_Holm[1:3]
priorct$model <- "M3_reference"; priorct$N <- nrow(d); priorct$P_Holm <- priorhm$P_Holm[4:5]
priorct$CI_lower <- priorct$beta_difference - qt(.975, priorct$df) * priorct$SE
priorct$CI_upper <- priorct$beta_difference + qt(.975, priorct$df) * priorct$SE
results$M3_reference <- list(coefficients = priorcf, contrasts = priorct, holm = priorhm)
ma <- fit("B2_2a_neutrophil_only", as.formula(paste("depression ~ log_Neutrophil +", covars)), ds, quasibinomial())
single <- extract(ma, "log_Neutrophil", "B2_2a_neutrophil_only")
csv(single, "B2_2a_neutrophil_only.csv")
run_joint("B2_2b_medications", ds, "+ any_nsaid + any_corticosteroid + any_statin")
domains <- list(c = subset(ds, !is.na(diabetes) & !is.na(hypertension)), d = subset(ds, !is.na(education) & !is.na(pir)), e = subset(ds, !is.na(alcohol_drinks)))
extras <- c(c = "+ diabetes + hypertension", d = "+ education + pir", e = "+ alcohol_drinks")
complete_counts <- data.frame()
for (letter in names(domains)) {
  sub <- domains[[letter]]; id <- paste0("B2_2", letter)
  run_joint(paste0(id, "_same_sample_M3"), sub)
  run_joint(paste0(id, "_extra_adjustment"), sub, extras[[letter]])
  complete_counts <- rbind(complete_counts, data.frame(task = id, N_total = nrow(d), N_complete = nrow(sub$variables), N_missing_any_additional = nrow(d) - nrow(sub$variables), missing_proportion = (nrow(d) - nrow(sub$variables)) / nrow(d)))
}
stopifnot(complete_counts$N_complete == c(33104, 29776, 22206))
run_joint("B2_2f_PHQ15", ds, outcome = "depression_strict")
stroke_ds <- subset(ds, stroke == 0); stopifnot(nrow(stroke_ds$variables) == 31245L)
run_joint("B2_2g_stroke0", stroke_ds)
run_joint("B2_2h_neutrophil_le9600", subset(ds, Neutrophil_count <= 9600))
run_joint("B2_2i_lymphocyte_le5000", subset(ds, Lymphocyte_count <= 5000))
csv(complete_counts, "B2_2_complete_case_counts.csv")
allcf <- do.call(rbind, lapply(results, `[[`, "coefficients")); allct <- do.call(rbind, lapply(results, `[[`, "contrasts")); allhm <- do.call(rbind, lapply(results, `[[`, "holm"))
csv(allcf, "B2_2_joint_coefficients.csv"); csv(allct, "B2_2_t_contrasts.csv"); csv(allhm, "B2_2_holm_five_tests.csv")
metrics <- c("N", unlist(lapply(cells, function(x) paste(x, c("OR", "CI_lower", "CI_upper", "P_raw", "P_Holm"), sep = "_"))), unlist(lapply(c("Neutrophil-Lymphocyte", "Neutrophil-Platelet"), function(x) paste(x, c("P_raw", "P_Holm", "df"), sep = "_"))))
side <- data.frame(metric = metrics)
for (id in names(results)) {
  r <- results[[id]]
  side[[id]] <- c(r$coefficients$N[1], unlist(lapply(cells, function(t) as.numeric(r$coefficients[r$coefficients$term == t, c("OR", "CI_lower", "CI_upper", "P_raw", "P_Holm")]))), unlist(lapply(c("Neutrophil-Lymphocyte", "Neutrophil-Platelet"), function(t) as.numeric(r$contrasts[r$contrasts$comparison == t, c("P_raw", "P_Holm", "df")]))))
}
side$B2_2a_neutrophil_only <- c(single$N, as.numeric(single[1, c("OR", "CI_lower", "CI_upper", "P_raw")]), NA, rep(NA_real_, length(metrics) - 6))
csv(side, "B2_2_side_by_side.csv")
for (letter in c("c", "d", "e")) {
  csv(side[, c("metric", paste0("B2_2", letter, "_same_sample_M3"), paste0("B2_2", letter, "_extra_adjustment"))], paste0("B2_2", letter, "_paired_same_sample.csv"))
}
ca <- so <- d
ca$outcome_type <- "PHQ9_CA"; ca$score <- ca$PHQ9_CA
so$outcome_type <- "PHQ9_SO"; so$score <- so$PHQ9_SO
long <- rbind(ca, so); long$outcome_type <- factor(long$outcome_type, levels = c("PHQ9_CA", "PHQ9_SO"))
dl <- svydesign(id = ~SDMVPSU + SEQN, strata = ~SDMVSTRA, weights = ~WTMEC_adj, nest = TRUE, data = long)
mh <- fit("B2_3_stacked_quasiPoisson_M3", score ~ outcome_type * (log_Neutrophil + age + sex + bmi + smoking + cycle_factor + race_eth), dl, quasipoisson(link = "log"))
it <- grep("outcome_typePHQ9_SO:log_Neutrophil|log_Neutrophil:outcome_typePHQ9_SO", names(coef(mh)), value = TRUE); stopifnot(length(it) == 1)
linear <- function(label, members) {
  v <- setNames(rep(0, length(coef(mh))), names(coef(mh))); v[members] <- 1
  obj <- svycontrast(mh, v); b <- as.numeric(coef(obj)); se <- as.numeric(SE(obj)); t <- b / se; df <- degf(dl)
  data.frame(result = label, N_persons = nrow(d), N_stacked = nobs(mh), beta_log_ratio = b, SE = se, ratio = exp(b), CI_lower = exp(b - qt(.975, df) * se), CI_upper = exp(b + qt(.975, df) * se), t = t, df = df, P_t = 2 * pt(-abs(t), df))
}
qp <- rbind(linear("SO_rate_ratio_per_ln", c("log_Neutrophil", it)), linear("CA_rate_ratio_per_ln", "log_Neutrophil"), linear("SO_to_CA_ratio_of_rate_ratios", it))
csv(qp, "B2_3_quasiPoisson_symptom_ratios.csv")
csv(registry[registry$model == "B2_3_stacked_quasiPoisson_M3", c("model", "converged", "iterations", "warning_count")], "B2_3_convergence.csv")
binary <- data.frame()
for (resp in c("SO_any", "CA_any")) {
  id <- paste0("B2_4_", resp)
  m <- fit(id, as.formula(paste(resp, "~ log_Neutrophil +", covars)), ds, quasibinomial())
  r <- extract(m, "log_Neutrophil", id); r$outcome_positive_N <- sum(d[[resp]] == 1)
  binary <- rbind(binary, r)
}
csv(binary, "B2_4_any_symptom_logistic.csv")
csv(registry, "B2_remaining_model_registry.csv"); csv(warnings_all, "B2_model_warnings.csv")
initial <- neutro_read_csv(file.path(out, "B2_initial_model_registry.csv"))
initial$link <- "logit"; initial$converged <- TRUE; initial$iterations <- NA_integer_; initial$warning_count <- NA_integer_
# Initial run predates explicit warning metadata; its original log is retained.
initial <- initial[, names(registry)]
csv(rbind(initial, registry), "B2_all_model_registry.csv")
csv(data.frame(status = "COMPLETE_WAITING_USER", remaining_models_run = nrow(registry), total_models_run = nrow(registry) + nrow(initial)), "B2_completion_status.csv")
neutro_cat("END COMPLETE_WAITING_USER REMAINING_MODELS_RUN=", nrow(registry), "\n")
