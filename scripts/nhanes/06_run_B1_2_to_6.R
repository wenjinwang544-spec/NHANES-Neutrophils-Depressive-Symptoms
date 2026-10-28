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
neutro_domain_entry('scripts/nhanes/06_run_B1_2_to_6.R')
# Derived candidate; see provenance/script_selection.csv. Invoke via run_job.py.
options(warn = 1)
suppressPackageStartupMessages(library(survey))
stopifnot(requireNamespace("Hmisc", quietly = TRUE))
out <- neutro_normalize("work/symptoms_base")
input <- neutro_normalize("work/cohort/supplement_input.rds")
csv <- function(x, name) {
  p <- file.path(out, name); if (neutro_file_exists(p)) stop("Output exists: ", p)
  neutro_write_csv(x, p, row.names = FALSE, na = "")
}
neutro_cat("START ", format(Sys.time(), tz = "Asia/Shanghai"), "\nINPUT ", input, "\nMD5 ", unname(neutro_md5(input)), "\n")
neutro_cat("已有 R 脚本只将编码 5 归入 Other，与 RDS 不一致，以 RDS 为准。race_eth 未重新编码。\n")
neutro_cat("对比检验类型=t，df=degf(design)，各设计分别计算；Holm 使用三系数P和两个对比P的完整精度。\n")
neutro_cat("Coefficient P/CI use model residual df. RCS exposure scale=log_Neutrophil, as existing Figure3 code; raw symptom outcomes in RCS.\n")
neutro_cat("Mutual adjustment: outcome z-standardized, other symptom predictor remains raw as specified.\n")
print(sessionInfo()); print(options()[c("contrasts", "survey.lonely.psu", "survey.adjust.domain.lonely")])
d <- neutro_readRDS(input)
required <- c("depression", "log_Neutrophil", "log_Lymphocyte", "log_Platelet", "age", "sex", "bmi", "smoking", "SDDSRVYR", "race_eth", "PHQ9_SO", "PHQ9_CA", "Neutrophil_count", "WTMEC_adj", "SDMVPSU", "SDMVSTRA", "SEQN")
stopifnot(all(required %in% names(d)), !anyNA(d[required]), nrow(d) == 33174L, sum(d$depression == 1) == 2902L, levels(d$smoking)[1] == "Never", !anyDuplicated(d$SEQN))
d$cycle_factor <- factor(d$SDDSRVYR, levels = 4:10)
race5 <- ifelse(d$RIDRETH3 %in% c(1, 2), "Hispanic", ifelse(d$RIDRETH3 == 3, "Non-Hispanic White", ifelse(d$RIDRETH3 == 4, "Non-Hispanic Black", ifelse(d$RIDRETH3 == 6, "Non-Hispanic Asian", ifelse(d$RIDRETH3 == 7, "Other/Multiracial", NA_character_)))))
d$race5 <- factor(race5, levels = c("Non-Hispanic White", "Hispanic", "Non-Hispanic Black", "Non-Hispanic Asian", "Other/Multiracial"))
ds <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC_adj, nest = TRUE, data = d)
registry <- data.frame()
fit <- function(id, formula, design, family) {
  neutro_cat("MODEL_START ", id, " N=", nrow(design$variables), " formula=", paste(deparse(formula), collapse = " "), " family=", family$family, "\n", sep = "")
  print(design$call)
  m <- svyglm(formula, design = design, family = family)
  print(summary(m)); print(confint(m))
  registry <<- rbind(registry, data.frame(model = id, formula = paste(deparse(formula), collapse = " "), N = nobs(m), family = family$family, df_residual = m$df.residual, df_design = degf(design)))
  neutro_cat("MODEL_END ", id, "\n")
  m
}
cells <- c("log_Neutrophil", "log_Lymphocyte", "log_Platelet")
base <- "age + sex + bmi + smoking"
adjustments <- c(M0 = base, M1 = paste(base, "+ cycle_factor"), M2 = paste(base, "+ race_eth"), M3 = paste(base, "+ cycle_factor + race_eth"))
coef_rows <- function(m, terms, id, exponentiate = TRUE) {
  b <- coef(m)[terms]; se <- sqrt(diag(vcov(m)))[terms]; ci <- confint(m)[terms, , drop = FALSE]; cs <- coef(summary(m))[terms, , drop = FALSE]
  x <- data.frame(model = id, N = nobs(m), term = terms, beta = unname(b), SE = unname(se), statistic = unname(b / se), df_residual = m$df.residual, P_raw = unname(cs[, ncol(cs)]), CI_lower_beta = unname(ci[, 1]), CI_upper_beta = unname(ci[, 2]), row.names = NULL)
  if (exponentiate) { x$OR <- exp(x$beta); x$CI_lower_OR <- exp(x$CI_lower_beta); x$CI_upper_OR <- exp(x$CI_upper_beta) }
  x
}
joint_result <- function(m, des, id) {
  cf <- coef_rows(m, cells, id)
  ctr <- function(expr, label) {
    o <- svycontrast(m, expr); b <- as.numeric(coef(o)); se <- as.numeric(SE(o)); t <- b / se; df <- degf(des)
    data.frame(model = id, N = nobs(m), comparison = label, beta_difference = b, SE = se, t = t, df = df, P_raw = 2 * pt(-abs(t), df), CI_lower = b - qt(.975, df) * se, CI_upper = b + qt(.975, df) * se)
  }
  ct <- rbind(ctr(quote(log_Neutrophil - log_Lymphocyte), "Neutrophil-Lymphocyte"), ctr(quote(log_Neutrophil - log_Platelet), "Neutrophil-Platelet"))
  hm <- data.frame(model = id, test = c(cells, ct$comparison), P_raw = c(cf$P_raw, ct$P_raw))
  hm$P_Holm <- p.adjust(hm$P_raw, method = "holm")
  cf$P_Holm <- hm$P_Holm[1:3]; ct$P_Holm <- hm$P_Holm[4:5]
  neutro_cat("SVYCONTRAST ", id, " t df=", degf(des), "\n"); print(ct); print(hm)
  list(coefficients = cf, contrasts = ct, holm = hm)
}
models <- results <- list()
for (id in names(adjustments)) {
  f <- as.formula(paste("depression ~", paste(cells, collapse = " + "), "+", adjustments[[id]]))
  models[[id]] <- fit(id, f, ds, quasibinomial())
  results[[id]] <- joint_result(models[[id]], ds, id)

}
allcf <- do.call(rbind, lapply(results, `[[`, "coefficients"))
allct <- do.call(rbind, lapply(results, `[[`, "contrasts"))
csv(allcf, "B1_2_joint_coefficients.csv"); csv(allct, "B1_2_t_contrasts.csv")
csv(do.call(rbind, lapply(results, `[[`, "holm")), "B1_2_holm_five_tests.csv")
wide <- function(res) {
  metrics <- c("N", unlist(lapply(cells, function(t) paste(t, c("OR", "CI_lower_OR", "CI_upper_OR", "P_raw", "P_Holm"), sep = "_"))), unlist(lapply(res[[1]]$contrasts$comparison, function(t) paste(t, c("beta_difference", "SE", "P_raw", "P_Holm", "df"), sep = "_"))))
  tb <- data.frame(metric = metrics)
  for (id in names(res)) {
    rr <- res[[id]]; vals <- c(rr$coefficients$N[1], unlist(lapply(cells, function(t) as.numeric(rr$coefficients[rr$coefficients$term == t, c("OR", "CI_lower_OR", "CI_upper_OR", "P_raw", "P_Holm")]))), unlist(lapply(rr$contrasts$comparison, function(t) as.numeric(rr$contrasts[rr$contrasts$comparison == t, c("beta_difference", "SE", "P_raw", "P_Holm", "df")]))))
    tb[[id]] <- vals
  }
  tb
}
csv(wide(results), "B1_2_M0_M3_side_by_side.csv")
logsd <- sqrt(as.numeric(coef(svyvar(~log_Neutrophil, ds))))
iqr <- as.numeric(coef(svyquantile(~Neutrophil_count, ds, c(.25, .75), ci = FALSE, qrule = "math")))
scale_multiplier <- c(1, logsd, log(5000 / 4000), log(iqr[2] / iqr[1]))
beta <- coef(models$M3)["log_Neutrophil"]; se <- sqrt(vcov(models$M3)["log_Neutrophil", "log_Neutrophil"]); df <- models$M3$df.residual
csv(data.frame(scale = c("per_1_ln_unit", "per_1_weighted_SD_log", "4000_to_5000_cells_uL", "weighted_Q25_to_Q75_cells_uL"), log_count_increment = scale_multiplier, count_from = c(NA, NA, 4000, iqr[1]), count_to = c(NA, NA, 5000, iqr[2]), beta_increment = as.numeric(beta * scale_multiplier), SE_increment = as.numeric(se * scale_multiplier), OR = as.numeric(exp(beta * scale_multiplier)), CI_lower = as.numeric(exp((beta - qt(.975, df) * se) * scale_multiplier)), CI_upper = as.numeric(exp((beta + qt(.975, df) * se) * scale_multiplier)), P_t = 2 * pt(-abs(beta / se), df), df_residual = df), "B1_3_M3_effect_scales.csv")
symmeans <- coef(svymean(~PHQ9_SO + PHQ9_CA, ds))
symsds <- sqrt(diag(as.matrix(svyvar(~PHQ9_SO + PHQ9_CA, ds))))
diagnostics <- do.call(rbind, lapply(c("PHQ9_SO", "PHQ9_CA"), function(v) {
  x <- d[[v]]; w <- d$WTMEC_adj / sum(d$WTMEC_adj); mu <- sum(w * x)
  q <- as.numeric(coef(svyquantile(as.formula(paste0("~", v)), ds, c(.25, .5, .75), ci = FALSE, qrule = "math")))
  m2 <- sum(w * (x - mu)^2); m3 <- sum(w * (x - mu)^3)
  data.frame(variable = v, N = length(x), weighted_mean = symmeans[v], weighted_SD = symsds[v], weighted_Q25 = q[1], weighted_median = q[2], weighted_Q75 = q[3], weighted_IQR = q[3] - q[1], weighted_zero_proportion = sum(w * (x == 0)), unweighted_zero_proportion = mean(x == 0), weighted_skewness_moment = m3 / m2^1.5, maximum = max(x), row.names = NULL)
}))
csv(diagnostics, "B1_4a_symptom_diagnostics.csv")
neutro_cat("DESCRIPTIVE skewness=sum(w*(x-mu)^3)/sum(w*(x-mu)^2)^1.5; weighted qrule=math; SD=svyvar.\n")
d$PHQ9_SO_z <- (d$PHQ9_SO - symmeans["PHQ9_SO"]) / symsds["PHQ9_SO"]
d$PHQ9_CA_z <- (d$PHQ9_CA - symmeans["PHQ9_CA"]) / symsds["PHQ9_CA"]
ds <- update(ds, PHQ9_SO_z = d$PHQ9_SO_z, PHQ9_CA_z = d$PHQ9_CA_z)
ca <- so <- d
ca$outcome_type <- "PHQ9_CA_z"; ca$score <- ca$PHQ9_CA_z
so$outcome_type <- "PHQ9_SO_z"; so$score <- so$PHQ9_SO_z
long <- rbind(ca, so); long$outcome_type <- factor(long$outcome_type, levels = c("PHQ9_CA_z", "PHQ9_SO_z"))
dl <- svydesign(id = ~SDMVPSU + SEQN, strata = ~SDMVSTRA, weights = ~WTMEC_adj, nest = TRUE, data = long)
mh <- fit("B1_4b_stacked_M3", score ~ outcome_type * (log_Neutrophil + age + sex + bmi + smoking + cycle_factor + race_eth), dl, gaussian())
int <- grep("outcome_typePHQ9_SO_z:log_Neutrophil|log_Neutrophil:outcome_typePHQ9_SO_z", names(coef(mh)), value = TRUE); stopifnot(length(int) == 1)
symcontrast <- function(label, members) {
  vec <- setNames(rep(0, length(coef(mh))), names(coef(mh))); vec[members] <- 1
  obj <- svycontrast(mh, vec); b <- as.numeric(coef(obj)); se <- as.numeric(SE(obj)); df <- degf(dl)
  data.frame(dimension = label, beta = b, SE = se, t = b / se, df = df, P_t_design_df = 2 * pt(-abs(b / se), df), CI_lower = b - qt(.975, df) * se, CI_upper = b + qt(.975, df) * se, df_residual_reference = mh$df.residual, P_t_residual_reference = 2 * pt(-abs(b / se), mh$df.residual))
}
csv(rbind(symcontrast("somatic", c("log_Neutrophil", int)), symcontrast("cognitive_affective", "log_Neutrophil"), symcontrast("SO_minus_CA", int)), "B1_4b_stacked_symptom_contrasts.csv")
mutual <- list()
for (v in c("SO", "CA")) {
  other <- if (v == "SO") "CA" else "SO"
  f <- as.formula(paste0("PHQ9_", v, "_z ~ log_Neutrophil + PHQ9_", other, " + ", adjustments[["M3"]]))
  mm <- fit(paste0("B1_4c_mutual_", v), f, ds, gaussian())
  mutual[[v]] <- coef_rows(mm, "log_Neutrophil", paste0("B1_4c_mutual_", v), FALSE)
}
csv(do.call(rbind, mutual), "B1_4c_mutually_adjusted_symptoms.csv")
# Figure3-compatible RCS on log-count scale: four full-sample weighted knots.
probs <- c(.05, .35, .65, .95)
knots <- as.numeric(coef(svyquantile(~log_Neutrophil, ds, probs, ci = FALSE, qrule = "math")))
csv(data.frame(probability = probs, knot_log_count = knots, knot_cells_uL = exp(knots)), "B1_5_spline_knots.csv")
basis <- Hmisc::rcspline.eval(d$log_Neutrophil, knots = knots, inclx = TRUE)
rcsds <- update(ds, rcs_linear = basis[, 1], rcs_nl1 = basis[, 2], rcs_nl2 = basis[, 3])
points <- c(2000, 3000, 4000, 5000, 7000, 9000)
bp <- Hmisc::rcspline.eval(log(points), knots = knots, inclx = TRUE)
br <- bp[points == 4000, ]; diff_basis <- sweep(bp, 2, br, "-")
spline_tests <- spline_points <- data.frame()
for (type in c("primary", "supplementary")) {
  extra <- if (type == "supplementary") "+ log_Lymphocyte + log_Platelet" else ""
  for (response in c("depression", "PHQ9_SO", "PHQ9_CA")) {
    id <- paste0("B1_5_", type, "_", response)
    f <- as.formula(paste(response, "~ rcs_linear + rcs_nl1 + rcs_nl2 +", adjustments[["M3"]], extra))
    fam <- if (response == "depression") quasibinomial() else gaussian()
    m <- fit(id, f, rcsds, fam)
    wt <- regTermTest(m, ~rcs_nl1 + rcs_nl2, method = "Wald", df = degf(rcsds))
    neutro_cat("NONLINEAR_WALD ", id, "\n"); print(wt)
    spline_tests <- rbind(spline_tests, data.frame(model = id, type = type, response = response, N = nobs(m), Wald_F = as.numeric(wt$Ftest), numerator_df = wt$df, denominator_df = degf(rcsds), P_nonlinear = as.numeric(wt$p)))
    terms <- c("rcs_linear", "rcs_nl1", "rcs_nl2"); b <- coef(m)[terms]; V <- vcov(m)[terms, terms]
    est <- as.numeric(diff_basis %*% b); se <- sqrt(pmax(0, rowSums((diff_basis %*% V) * diff_basis)))
    df <- degf(rcsds); lo <- est - qt(.975, df) * se; hi <- est + qt(.975, df) * se
    rows <- data.frame(model = id, type = type, response = response, count_cells_uL = points, reference_cells_uL = 4000, estimate_link_difference = est, SE_link_difference = se, df = df, scale = if (response == "depression") "OR_vs_4000" else "raw_symptom_score_difference_vs_4000", estimate = if (response == "depression") exp(est) else est, CI_lower = if (response == "depression") exp(lo) else lo, CI_upper = if (response == "depression") exp(hi) else hi)
    spline_points <- rbind(spline_points, rows)
  }
}
csv(spline_tests, "B1_5_nonlinearity_Wald.csv"); csv(spline_points, "B1_5_point_estimates.csv")
low <- do.call(rbind, lapply(levels(d$race_eth), function(r) {
  ix <- d$race_eth == r; small <- ix & d$Neutrophil_count <= 3000
  sub <- subset(ds, race_eth == r & Neutrophil_count <= 3000)
  pr <- svymean(~depression, sub)
  data.frame(race_eth = r, N_race = sum(ix), N_le_3000 = sum(small), depression_N_le_3000 = sum(d$depression[small] == 1), depression_proportion_unweighted = mean(d$depression[small]), depression_proportion_weighted = as.numeric(coef(pr)), weighted_proportion_SE = as.numeric(SE(pr)))
}))
csv(low, "B1_5_low_count_by_race.csv")
# Explicit domain subset, never filter rows then rebuild the survey design.
subds <- subset(ds, SDDSRVYR >= 7)
stopifnot(!anyNA(subds$variables$race5), !anyNA(subds$variables$RIDRETH3))
csv(as.data.frame(table(race5 = subds$variables$race5, useNA = "always")), "B1_6_race5_counts.csv")
csv(data.frame(N = nrow(subds$variables), depression_N = sum(subds$variables$depression == 1), df_design = degf(subds), minimum_cycle = min(subds$variables$SDDSRVYR), maximum_cycle = max(subds$variables$SDDSRVYR)), "B1_6_subsample_counts.csv")
subresults <- list()
for (race in c("race5", "race_eth")) {
  id <- paste0("B1_6_", race)
  f <- as.formula(paste("depression ~", paste(cells, collapse = " + "), "+", base, "+ cycle_factor +", race))
  m <- fit(id, f, subds, quasibinomial())
  subresults[[race]] <- joint_result(m, subds, id)
}
csv(do.call(rbind, lapply(subresults, `[[`, "coefficients")), "B1_6_joint_coefficients.csv")
csv(do.call(rbind, lapply(subresults, `[[`, "contrasts")), "B1_6_t_contrasts.csv")
csv(wide(subresults), "B1_6_four_vs_five_race_side_by_side.csv")
csv(registry, "B1_2_to_6_model_registry.csv")
csv(data.frame(stage = "B1-2_to_6", status = "COMPLETE_WAITING_USER", models_run = nrow(registry)), "B1_completion_status.csv")
neutro_cat("END COMPLETE_WAITING_USER MODELS_RUN=", nrow(registry), "\n")
