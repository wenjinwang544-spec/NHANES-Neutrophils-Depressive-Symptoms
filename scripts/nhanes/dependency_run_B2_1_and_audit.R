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
neutro_domain_entry('scripts/nhanes/dependency_run_B2_1_and_audit.R')
options(warn = 1)
suppressPackageStartupMessages(library(survey))
out <- neutro_normalize("work/symptoms_extended")
input <- neutro_normalize("work/cohort/supplement_input.rds")
csv <- function(x, name) {
  p <- file.path(out, name); stopifnot(!neutro_file_exists(p))
  neutro_write_csv(x, p, row.names = FALSE, na = "")
}
neutro_cat("START ", format(Sys.time(), tz = "Asia/Shanghai"), "\nINPUT ", input, "\nMD5 ", unname(neutro_md5(input)), "\n")
print(sessionInfo())
neutro_cat("race_eth=RDS stored factor; cycle_factor=factor(SDDSRVYR); contrasts=t,df=degf(design).\n")
d <- neutro_readRDS(input)
d$cycle_factor <- factor(d$SDDSRVYR, levels = 4:10)
vs <- c("any_nsaid", "any_corticosteroid", "any_statin", "diabetes", "hypertension", "education", "pir", "alcohol_drinks", "stroke")
csv(data.frame(variable = vs, class = sapply(d[vs], function(x) paste(class(x), collapse = ";")), N_total = nrow(d), N_missing = sapply(d[vs], function(x) sum(is.na(x))), missing_proportion = sapply(d[vs], function(x) mean(is.na(x))), levels = sapply(d[vs], function(x) paste(levels(x), collapse = ";"))), "B2_variable_missing_audit.csv")
csv(data.frame(set = c("diabetes+hypertension", "education+pir", "alcohol_drinks"), N_complete = c(sum(complete.cases(d[, c("diabetes", "hypertension")])), sum(complete.cases(d[, c("education", "pir")])), sum(!is.na(d$alcohol_drinks))), N_total = nrow(d)), "B2_adjustment_set_complete_counts.csv")
csv(as.data.frame(table(stroke = d$stroke, useNA = "always")), "B2_stroke_counts.csv")
ds <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC_adj, nest = TRUE, data = d)
formula <- depression ~ log_Neutrophil + log_Lymphocyte + log_Platelet + age + sex + bmi + smoking + cycle_factor + race_eth
stopifnot(nrow(d) == 33174L, sum(d$depression == 1) == 2902L, !anyNA(model.frame(formula, d)))
neutro_cat("MODEL_START B2_1_M3 N=", nrow(d), " formula=", paste(deparse(formula), collapse = " "), " family=quasibinomial\n")
print(ds$call)
m <- svyglm(formula, ds, family = quasibinomial())
print(summary(m)); print(confint(m)); neutro_cat("MODEL_END B2_1_M3\n")
terms <- c("log_Neutrophil", "log_Lymphocyte", "log_Platelet")
ct <- coef(summary(m))[terms, , drop = FALSE]; ci <- confint(m)[terms, , drop = FALSE]
rows <- data.frame(term = terms, N = nobs(m), beta = unname(coef(m)[terms]), SE = unname(sqrt(diag(vcov(m)))[terms]), OR = unname(exp(coef(m)[terms])), CI_lower = unname(exp(ci[, 1])), CI_upper = unname(exp(ci[, 2])), P_raw = unname(ct[, ncol(ct)]), df_residual = m$df.residual)
contrast <- function(expr, label) {
  o <- svycontrast(m, expr); b <- as.numeric(coef(o)); se <- as.numeric(SE(o)); t <- b / se; df <- degf(ds)
  data.frame(comparison = label, beta_difference = b, SE = se, t = t, df = df, P_raw = 2 * pt(-abs(t), df))
}
contrasts <- rbind(contrast(quote(log_Neutrophil - log_Lymphocyte), "Neutrophil-Lymphocyte"), contrast(quote(log_Neutrophil - log_Platelet), "Neutrophil-Platelet"))
holm <- data.frame(test = c(terms, contrasts$comparison), P_raw = c(rows$P_raw, contrasts$P_raw))
holm$P_Holm <- p.adjust(holm$P_raw, "holm")
neutro_cat("SVYCONTRAST t df=", degf(ds), "\n"); print(contrasts); print(holm)
csv(rows, "B2_1_M3_joint_coefficients.csv"); csv(contrasts, "B2_1_M3_t_contrasts.csv"); csv(holm, "B2_1_M3_holm_five_tests.csv")
old <- neutro_read_csv("work/symptoms_base/B1_2_joint_coefficients.csv"); old <- old[old$model == "M3", ]
oldh <- neutro_read_csv("work/symptoms_base/B1_2_holm_five_tests.csv"); oldh <- oldh[oldh$model == "M3", ]
validation <- rbind(data.frame(item = terms, metric = "beta", previous = old$beta, current = rows$beta), data.frame(item = holm$test, metric = "P_Holm", previous = oldh$P_Holm, current = holm$P_Holm))
validation$difference <- validation$current - validation$previous; validation$pass <- abs(validation$difference) <= 1e-12
csv(validation, "B2_1_M3_validation.csv")
csv(data.frame(model = "B2_1_M3", formula = paste(deparse(formula), collapse = " "), N = nobs(m), family = "quasibinomial", df_residual = m$df.residual, df_design = degf(ds)), "B2_initial_model_registry.csv")
stopifnot(all(validation$pass))
csv(data.frame(stage = "B2-1", status = "WAITING_SPECIFICATION_CONFIRMATION", models_run = 1), "B2_initial_status.csv")
neutro_cat("END WAITING_SPECIFICATION_CONFIRMATION MODELS_RUN=1\n")
