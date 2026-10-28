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
neutro_domain_entry('scripts/mr/upstream_online_MR_UNVALIDATED.R')
input_root<-getwd()
# 验证是否有效
library(ieugwasr)
# User account lookup removed from candidate code.

# ============================================================================
# 研究题目：Causal Link Between Neutrophil-Driven Systemic Inflammation
#           and Depression Risk
# 结局数据：PGC MDD2025 EUR no23andMe-noUKBB                    ← ✅ 已修正
#           (pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11_formatted.tsv.gz)
# 修正原因：暴露GWAS(ieu-b-34/ieu-b-32)源自FinnGen+UKBB，
#           原版no23andMe结局GWAS含UKBB，存在样本重叠，违反双样本MR独立性假设
#           noUKBB版本明确剔除UK Biobank，消除样本重叠偏倚
# 暴露数据：中性粒细胞(ieu-b-34) / 淋巴细胞(ieu-b-32) / 血小板(ukb-d-30080_irnt)
# ============================================================================

# 设置输出目录（所有CSV/PNG文件都会保存到这里）
output_dir <- neutro_normalize('work/mr_online'); neutro_setwd(output_dir)  # 改成你想要的路径

neutro_cat("=== 完整MR分析脚本（已修复列名） ===\n")
neutro_cat("基于本地PGC MDD2025数据\n")
neutro_cat("文件：pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11_formatted.tsv.gz\n\n")

# ============================================================================
# 第1步：加载必要的包
# ============================================================================
# install.packages(c("TwoSampleMR","data.table","dplyr","ggplot2","ieugwasr","MVMR"))
neutro_cat("【第1步】加载R包...\n")

library(TwoSampleMR)
library(data.table)
library(dplyr)
library(ggplot2)
library(MVMR)

neutro_cat("✓ 包加载成功\n\n")

# ============================================================================
# 第二步：读取PGC MDD2025结局数据
# ============================================================================



# ============================================================================
# 第2步：读取PGC MDD2025结局数据
# ============================================================================

neutro_cat("【第2步】读取PGC MDD2025数据（noUKBB版本）...\n\n")

# ✅ 已切换至 noUKBB 版本，解决暴露GWAS(UKBB)与结局GWAS样本重叠问题
# 原文件：pgc-mdd2025_no23andMe_eur_v3-49-24-11_formatted.tsv.gz        ← 含UKBB ⚠️
# 新文件：pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11_formatted.tsv.gz ← 不含UKBB ✅
pgc_file <- file.path(input_root, 'data/gwas/pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11_formatted.tsv.gz')
# 或直接指定路径：
# pgc_file <- "/你的路径/pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11_formatted.tsv.gz"

if (!neutro_file_exists(pgc_file)) stop("✗ 找不到文件，请检查路径")

neutro_cat("正在读取数据（约需1-5分钟）...\n")
mdd_data <- neutro_fread(pgc_file, showProgress = TRUE)
neutro_cat("✓ 读取完成：", nrow(mdd_data), "行 /", ncol(mdd_data), "列\n\n")

# ── 列名映射（daner格式：ID/EA/NEA/BETA/SE/PVAL/FCAS/NCAS/NCON/NEFF）────────
snp_col      <- "ID"
pval_col     <- "PVAL"
beta_col     <- "BETA"
se_col       <- "SE"
ea_col       <- "EA"
nea_col      <- "NEA"
eaf_col      <- "FCAS"      # 病例组等位基因频率（也可用FCON）
ncase_col    <- "NCAS"
ncontrol_col <- "NCON"
n_col        <- "NEFF"      # 有效样本量

neutro_cat("列名映射确认：\n")
neutro_cat("  SNP:      ", snp_col,      "\n")
neutro_cat("  BETA:     ", beta_col,     "\n")
neutro_cat("  SE:       ", se_col,       "\n")
neutro_cat("  P值:      ", pval_col,     "\n")
neutro_cat("  效应等位: ", ea_col,       "\n")
neutro_cat("  另一等位: ", nea_col,      "\n")
neutro_cat("  EAF:      ", eaf_col,      "\n")
neutro_cat("  病例数:   ", ncase_col,    "\n")
neutro_cat("  对照数:   ", ncontrol_col, "\n")
neutro_cat("  有效N:    ", n_col,        "\n\n")

# ============================================================================
# 第3步：清理并格式化结局数据
# ============================================================================

neutro_cat("【第3步】清理并格式化结局数据...\n\n")

# 使用自动检测到的列名过滤缺失值
mdd_clean <- mdd_data[
  !is.na(mdd_data[[snp_col]])    &
    !is.na(mdd_data[[beta_col]]) &
    !is.na(mdd_data[[se_col]])   &
    !is.na(mdd_data[[ea_col]])   &
    !is.na(mdd_data[[nea_col]])  &
    !is.na(mdd_data[[pval_col]]), ]

neutro_cat("  原始行数：", nrow(mdd_data), "\n")
neutro_cat("  清理后行数：", nrow(mdd_clean), "\n")
neutro_cat("  移除行数：", nrow(mdd_data) - nrow(mdd_clean), "\n\n")

if (nrow(mdd_clean) == 0) {
  neutro_cat("✗ 清理后0行，诊断信息：\n")
  neutro_cat("  SNP列(", snp_col, ")前5值：", head(mdd_data[[snp_col]], 5), "\n")
  neutro_cat("  beta前5值：",            head(mdd_data[["beta"]], 5), "\n")
  stop("数据清理失败，请检查上方列名和前5行内容")
}

# 构建格式化数据框（使用检测到的列名）
format_df <- data.frame(
  SNP  = mdd_clean[[snp_col]],
  BETA = as.numeric(mdd_clean[[beta_col]]),
  SE   = as.numeric(mdd_clean[[se_col]]),
  A1   = mdd_clean[[ea_col]],
  A2   = mdd_clean[[nea_col]],
  P    = as.numeric(mdd_clean[[pval_col]]),
  FRQ        = as.numeric(mdd_clean[[eaf_col]]),
  N_CASES    = as.numeric(mdd_clean[[ncase_col]]),
  N_CONTROLS = as.numeric(mdd_clean[[ncontrol_col]]),
  N_TOTAL    = as.numeric(mdd_clean[[n_col]])
)

# 仅对核心列过滤缺失值
format_df <- format_df[
  !is.na(format_df$SNP)  &
    !is.na(format_df$BETA) &
    !is.na(format_df$SE)   &
    !is.na(format_df$A1)   &
    !is.na(format_df$A2)   &
    !is.na(format_df$P), ]

neutro_cat("  格式化后行数：", nrow(format_df), "\n\n")

# 格式化为TwoSampleMR标准结局数据（按实际存在的列动态传参）
format_args <- list(
  dat               = format_df,
  type              = "outcome",
  snp_col           = "SNP",
  beta_col          = "BETA",
  se_col            = "SE",
  effect_allele_col = "A1",
  other_allele_col  = "A2",
  pval_col          = "P"
)
if ("FRQ"        %in% names(format_df)) format_args$eaf_col        <- "FRQ"
if ("N_CASES"    %in% names(format_df)) format_args$ncase_col      <- "N_CASES"
if ("N_CONTROLS" %in% names(format_df)) format_args$ncontrol_col   <- "N_CONTROLS"
if ("N_TOTAL"    %in% names(format_df)) format_args$samplesize_col <- "N_TOTAL"

outcome_data <- do.call(format_data, format_args)

# 设置结局名称
outcome_data$outcome <- "Major Depressive Disorder (PGC MDD2025 EUR)"

# 验证
neutro_cat("验证结局数据：\n")
neutro_cat("  SNP总数：      ", nrow(outcome_data), "\n")
neutro_cat("  结局名称：     ", unique(outcome_data$outcome), "\n")
neutro_cat("  含样本量列：   ", "samplesize.outcome" %in% colnames(outcome_data), "\n")
neutro_cat("  含病例数列：   ", "ncase.outcome"      %in% colnames(outcome_data), "\n")
neutro_cat("  平均有效N：    ", round(mean(outcome_data$samplesize.outcome, na.rm=TRUE)), "\n")
neutro_cat("  平均病例数：   ", round(mean(outcome_data$ncase.outcome,      na.rm=TRUE)), "\n\n")

# 释放内存
rm(mdd_data, mdd_clean, format_df)
gc()

# ============================================================================
# 第4步：从IEU OpenGWAS提取暴露数据（工具变量）
# ============================================================================
library(TwoSampleMR)
library(ieugwasr)

# 设置Token
options(ieugwasr_token = Sys.getenv("IEUGWAS_TOKEN"))

# 定义暴露变量（同源Vuckovic 2020数据）
exposures <- list(
  neutrophil = "ieu-b-34",         # 中性粒细胞
  lymphocyte = "ieu-b-32",         # 淋巴细胞
  platelet   = "ukb-d-30080_irnt"  # 血小板
)

exposure_data_list <- list()

for (exp_name in names(exposures)) {
  neutro_cat("提取", exp_name, "数据...\n")
  
  tryCatch({
    exp_dat <- extract_instruments(
      outcomes = exposures[[exp_name]],
      p1       = 5e-8,    # 全基因组显著性阈值
      clump    = TRUE,    # LD剪枝
      r2       = 0.001,   # LD r²阈值（严格）
      kb       = 10000    # 窗口10Mb
    )
    
    if (!is.null(exp_dat) && nrow(exp_dat) > 0) {
      
      # ✅ 补充：过滤弱工具变量（F统计量 < 10）
      exp_dat$F_stat <- (exp_dat$beta.exposure / exp_dat$se.exposure)^2
      weak_iv <- sum(exp_dat$F_stat < 10)
      if (weak_iv > 0) {
        neutro_cat("  ⚠ 移除", weak_iv, "个弱工具变量(F<10)\n")
        exp_dat <- exp_dat[exp_dat$F_stat >= 10, ]
      }
      
      exp_dat$exposure <- exp_name
      exposure_data_list[[exp_name]] <- exp_dat
      
      neutro_cat("  ✓ 成功提取", nrow(exp_dat), "个SNP\n")
      neutro_cat("  ✓ 平均F统计量:", round(mean(exp_dat$F_stat), 1), "\n")
      
    } else {
      neutro_cat("  ✗ 无法提取数据（返回为空）\n")
    }
    
  }, error = function(e) {
    neutro_cat("  ✗ 错误:", e$message, "\n")
  })
  
  neutro_cat("\n")
}

# 检查结果
if (length(exposure_data_list) == 0) {
  stop("✗ 错误: 无法提取任何暴露数据，请检查Token和GWAS ID")
}

neutro_cat("✓ 暴露数据提取完成\n")
neutro_cat("  成功提取的变量:", paste(names(exposure_data_list), collapse = ", "), "\n")
neutro_cat("  各变量SNP数:\n")
for (n in names(exposure_data_list)) {
  neutro_cat("    -", n, ":", nrow(exposure_data_list[[n]]), "个SNP\n")
}

# 在保存前添加去重步骤
for (n in names(exposure_data_list)) {
  dat <- exposure_data_list[[n]]
  dup_snps <- dat$SNP[duplicated(dat$SNP)]
  if (length(dup_snps) > 0) {
    neutro_cat("去除", n, "中的重复SNP:", paste(dup_snps, collapse=", "), "\n")
    exposure_data_list[[n]] <- dat[!duplicated(dat$SNP), ]
  }
}

# ✅ 立即保存，避免重复提取
neutro_saveRDS(exposure_data_list, "exposure_data_list.rds")
neutro_cat("\n✓ 已保存至 exposure_data_list.rds\n")

neutro_cat("✓ 暴露数据提取完成\n")
neutro_cat("  成功提取的变量:", paste(names(exposure_data_list), collapse = ", "), "\n\n")

# ============================================================================
# 第5步：检查暴露数据 + 执行数据调和
# ============================================================================

neutro_cat("【第5步】检查暴露数据...\n\n")

if (!exists("exposure_data_list")) {
  neutro_cat("✗ 错误：找不到 'exposure_data_list' 变量\n\n")
  neutro_cat("选项1：从保存的文件加载\n")
  neutro_cat("  exposure_data_list <- readRDS('exposure_data_list.rds')\n\n")
  neutro_cat("选项2：重新提取暴露数据\n")
  stop("请先加载或提取暴露数据")
}

neutro_cat("✓ 已检测到暴露数据\n")
neutro_cat("  暴露变量数：", length(exposure_data_list), "\n")
neutro_cat("  变量名：", paste(names(exposure_data_list), collapse = ", "), "\n\n")

# --- 新增：去除重复SNP ---
for (n in names(exposure_data_list)) {
  dat <- exposure_data_list[[n]]
  dup_snps <- dat$SNP[duplicated(dat$SNP)]
  if (length(dup_snps) > 0) {
    neutro_cat("⚠ 去除", n, "中的重复SNP:", paste(dup_snps, collapse=", "), "\n")
    exposure_data_list[[n]] <- dat[!duplicated(dat$SNP), ]
  }
}

# --- 新增：执行数据调和 ---
neutro_cat("开始数据调和...\n\n")

harmonised_list <- list()

for (exp_name in names(exposure_data_list)) {
  neutro_cat("调和", exp_name, "数据...\n")
  
  tryCatch({
    harmonised <- harmonise_data(
      exposure_dat = exposure_data_list[[exp_name]],
      outcome_dat  = outcome_data,
      action = 2
    )
    
    # 仅保留可用SNP
    harmonised_keep <- harmonised[harmonised$mr_keep == TRUE, ]
    
    neutro_cat("  原始SNP数：", nrow(harmonised), "\n")
    neutro_cat("  可用SNP数：", nrow(harmonised_keep), "\n")
    neutro_cat("  剔除SNP数：", nrow(harmonised) - nrow(harmonised_keep), 
        "（回文/不匹配）\n\n")
    
    harmonised_list[[exp_name]] <- harmonised_keep
    
  }, error = function(e) {
    neutro_cat("  ✗ 调和失败：", e$message, "\n\n")
  })
}

# 检查调和结果
if (length(harmonised_list) == 0) {
  stop("✗ 所有暴露变量调和失败")
}
neutro_saveRDS(harmonised_list, "harmonised_list.rds")
neutro_cat("✓ 数据调和完成，可进行MR分析\n\n")

# SNP概览（调和前后对比）

neutro_cat("── 调和前后SNP数量对比 ──\n")
for (exp_name in names(harmonised_list)) {
  n_before <- nrow(exposure_data_list[[exp_name]])
  n_after  <- nrow(harmonised_list[[exp_name]])
  neutro_cat(sprintf("  %-12s 提取:%d → 可用:%d（剔除%d个）\n",
              exp_name, n_before, n_after, n_before - n_after))
}
neutro_cat("\n")

# ============================================================================
# 第6步：UVMR主分析
# ============================================================================

mr_methods <- c(
  "mr_ivw",
  "mr_egger_regression",
  "mr_weighted_median",
  "mr_weighted_mode"
)

uvmr_results_list  <- list()
heterogeneity_list <- list()
pleiotropy_list    <- list()
steiger_list       <- list()

for (exp_name in names(harmonised_list)) {
  
  neutro_cat("════════════════════════════════════\n")
  neutro_cat("  暴露变量：", exp_name, "\n")
  neutro_cat("════════════════════════════════════\n")
  
  dat <- harmonised_list[[exp_name]]
  
  if (nrow(dat) < 3) {
    neutro_cat("  ✗ 可用SNP不足3个，跳过\n\n")
    next
  }
  
  # [1] 主MR分析（四种方法）
  neutro_cat("  [1] 主MR分析...\n")
  mr_res        <- mr(dat, method_list = mr_methods)
  mr_res$OR     <- exp(mr_res$b)
  mr_res$OR_lci <- exp(mr_res$b - 1.96 * mr_res$se)
  mr_res$OR_uci <- exp(mr_res$b + 1.96 * mr_res$se)
  uvmr_results_list[[exp_name]] <- mr_res
  print(mr_res[, c("method","nsnp","b","se","pval","OR","OR_lci","OR_uci")])
  neutro_cat("\n")
  
  # [2] 异质性检验（Cochran Q）
  neutro_cat("  [2] 异质性检验（Cochran Q）...\n")
  het_res <- mr_heterogeneity(dat,
                              method_list = c("mr_egger_regression", "mr_ivw"))
  heterogeneity_list[[exp_name]] <- het_res
  print(het_res[, c("method","Q","Q_df","Q_pval")])
  ivw_qp <- het_res$Q_pval[grepl("Inverse", het_res$method)]
  neutro_cat("  →", ifelse(ivw_qp < 0.05,
                    paste0("⚠ 存在显著异质性 (p=", formatC(ivw_qp, format="e", digits=2), ")"),
                    paste0("✓ 无显著异质性 (p=",   round(ivw_qp, 4), ")")), "\n\n")
  
  # [3] 水平多效性检验（MR-Egger截距）
  neutro_cat("  [3] MR-Egger多效性检验...\n")
  plei_res <- mr_pleiotropy_test(dat)
  pleiotropy_list[[exp_name]] <- plei_res
  neutro_cat("  → Egger截距：", round(plei_res$egger_intercept, 6),
      "  p =", round(plei_res$pval, 4),
      ifelse(plei_res$pval < 0.05,
             "  ⚠ 存在水平多效性",
             "  ✓ 无水平多效性"), "\n\n")
  
  # [4] Steiger方向检验（需样本量信息）
  neutro_cat("  [4] Steiger方向检验...\n")
  tryCatch({
    # 为二分类结局（MDD）计算准确的r值
    dat$r.outcome <- get_r_from_lor(
      lor        = dat$beta.outcome,
      af         = dat$eaf.outcome,
      ncase      = dat$ncase.outcome,
      ncontrol   = dat$ncontrol.outcome,
      prevalence = 0.15    # MDD终身患病率约15%
    )
    st_res <- directionality_test(dat)
    steiger_list[[exp_name]] <- st_res
    neutro_cat("  → 方向正确：", st_res$correct_causal_direction,
        "  p =", formatC(st_res$steiger_pval, format="e", digits=2),
        ifelse(!st_res$correct_causal_direction,
               "  ⚠ 注意反向因果风险", "  ✓"), "\n\n")
  }, error = function(e) {
    neutro_cat("  ⚠ Steiger检验跳过（缺少样本量信息）：", e$message, "\n\n")
  })
  
  # [5] 留一法稳健性检验
  neutro_cat("  [5] 留一法稳健性检验...\n")
  loo_res   <- mr_leaveoneout(dat)
  ivw_b     <- mr_res$b[grepl("Inverse", mr_res$method)]
  loo_range <- range(loo_res$b)
  neutro_cat("  → IVW效应值：  ", round(ivw_b, 4), "\n")
  neutro_cat("  → 留一法范围：[", round(loo_range[1], 4), ",",
      round(loo_range[2], 4), "]\n")
  neutro_cat("  →", ifelse(diff(loo_range) > 0.5 * abs(ivw_b),
                    "⚠ 存在有影响力SNP，建议查看留一法图",
                    "✓ 结果稳健，无单一SNP主导效应"), "\n\n")
}

# ============================================================================
# 第7步：UVMR结果汇总与保存
# ============================================================================

neutro_cat("【第7步】UVMR结果汇总...\n\n")

all_mr_res        <- do.call(rbind, uvmr_results_list)
all_mr_res$OR     <- exp(all_mr_res$b)
all_mr_res$OR_lci <- exp(all_mr_res$b - 1.96 * all_mr_res$se)
all_mr_res$OR_uci <- exp(all_mr_res$b + 1.96 * all_mr_res$se)

result_table <- data.frame(
  Exposure = all_mr_res$exposure,
  Method   = all_mr_res$method,
  nSNP     = all_mr_res$nsnp,
  Beta     = round(all_mr_res$b,   4),
  SE       = round(all_mr_res$se,  4),
  Pval     = formatC(all_mr_res$pval, format = "e", digits = 2),
  OR       = round(all_mr_res$OR,     3),
  CI_95    = paste0(round(all_mr_res$OR_lci, 3),
                    " – ",
                    round(all_mr_res$OR_uci, 3))
)

neutro_cat("── IVW主要结果 ──\n")
print(result_table[grepl("Inverse variance weighted", result_table$Method), ])

neutro_write_csv(result_table, "uvmr_results_full.csv", row.names = FALSE)
neutro_cat("\n✓ 完整UVMR结果已保存至 uvmr_results_full.csv\n\n")

# ============================================================================
# 第8步：可视化（散点图 / 留一法图 / 漏斗图）
# ============================================================================

neutro_cat("【第8步】生成可视化图表...\n\n")

for (exp_name in names(harmonised_list)) {
  
  neutro_cat("生成", exp_name, "图表...\n")
  dat <- harmonised_list[[exp_name]]
  
  # 过滤不完整SNP避免ggplot2警告
  dat_plot <- dat[
    !is.na(dat$beta.exposure) & !is.na(dat$beta.outcome) &
      !is.na(dat$se.exposure)   & !is.na(dat$se.outcome), ]
  
  tryCatch({
    
    # 散点图
    p_scatter <- mr_scatter_plot(
      mr(dat_plot, method_list = mr_methods), dat_plot)
    neutro_ggsave(paste0("scatter_", exp_name, ".png"),
           p_scatter[[1]], width = 7, height = 7, dpi = 300)
    
    # 留一法图
    p_loo <- mr_leaveoneout_plot(mr_leaveoneout(dat_plot))
    neutro_ggsave(paste0("loo_", exp_name, ".png"),
           p_loo[[1]], width = 7, height = 10, dpi = 300)
    
    # 漏斗图
    p_funnel <- mr_funnel_plot(mr_singlesnp(dat_plot))
    neutro_ggsave(paste0("funnel_", exp_name, ".png"),
           p_funnel[[1]], width = 7, height = 7, dpi = 300)
    
    neutro_cat("  ✓ 三张图已保存\n")
    
  }, error = function(e) {
    neutro_cat("  ✗ 绘图失败：", e$message, "\n")
  })
}
neutro_cat("\n")

# ============================================================================
# 第9步：MVMR分析（多变量孟德尔随机化）
# ============================================================================

neutro_cat("【第9步】MVMR分析...\n\n")

# --- 9.1 提取MVMR暴露数据 ---
neutro_cat("  [9.1] 提取MVMR暴露数据（合并三个GWAS的SNP）...\n")

mvmr_exposure_dat <- mv_extract_exposures(
  id             = c("ieu-b-34", "ieu-b-32", "ukb-d-30080_irnt"),
  clump_r2       = 0.001,
  clump_kb       = 10000,
  pval_threshold = 5e-8
)

neutro_cat("  ✓ MVMR暴露SNP总数：", nrow(mvmr_exposure_dat), "\n\n")

# --- 9.2 匹配结局数据 ---
neutro_cat("  [9.2] 匹配结局数据...\n")
mvmr_outcome_dat <- outcome_data[
  outcome_data$SNP %in% mvmr_exposure_dat$SNP, ]
neutro_cat("  ✓ 结局匹配SNP数：", nrow(mvmr_outcome_dat), "\n\n")

# --- 9.3 MVMR数据调和 ---
neutro_cat("  [9.3] MVMR数据调和...\n")
mvmr_harmonised <- mv_harmonise_data(
  exposure_dat = mvmr_exposure_dat,
  outcome_dat  = mvmr_outcome_dat
)
neutro_cat("  ✓ 调和后SNP数：", nrow(mvmr_harmonised$exposure_beta), "\n\n")

# --- 9.4 MVMR主分析（IVW）---
neutro_cat("  [9.4] MVMR主分析...\n")
mvmr_res              <- mv_multiple(mvmr_harmonised)
mvmr_res$result$OR    <- exp(mvmr_res$result$b)
mvmr_res$result$OR_lci <- exp(mvmr_res$result$b - 1.96 * mvmr_res$result$se)
mvmr_res$result$OR_uci <- exp(mvmr_res$result$b + 1.96 * mvmr_res$result$se)

neutro_cat("  MVMR结果（控制三个暴露后的独立效应）：\n")
print(mvmr_res$result[, c("exposure","nsnp","b","se","pval","OR","OR_lci","OR_uci")])
neutro_cat("\n")

# --- 9.5 MVMR敏感性分析 ---
neutro_cat("  [9.5] MVMR敏感性分析...\n\n")

mvmr_strength <- NULL

tryCatch({
  
  mvmr_format <- format_mvmr(
    BXGs   = mvmr_harmonised$exposure_beta,
    BYG    = mvmr_harmonised$outcome_beta,
    seBXGs = mvmr_harmonised$exposure_se,
    seBYG  = mvmr_harmonised$outcome_se,
    RSID   = rownames(mvmr_harmonised$exposure_beta)
  )
  
  # 条件F统计量（工具变量强度）
  neutro_cat("  条件F统计量（>10为强工具变量）：\n")
  mvmr_strength <- strength_mvmr(r_input = mvmr_format, gencov = 0)
  print(mvmr_strength)
  neutro_cat("\n")
  
  # Q统计量（多效性检验）
  neutro_cat("  Q统计量多效性检验（p>0.05为无显著多效性）：\n")
  pleio_res <- pleiotropy_mvmr(r_input = mvmr_format, gencov = 0)
  print(pleio_res)
  neutro_cat("\n")
  
}, error = function(e) {
  neutro_cat("  ⚠ MVMR敏感性分析失败：", e$message, "\n\n")
})

# 保存MVMR结果
neutro_write_csv(mvmr_res$result, "mvmr_results.csv", row.names = FALSE)
neutro_cat("  ✓ MVMR结果已保存至 mvmr_results.csv\n\n")

# ============================================================================
# 第10步：UVMR vs MVMR结果对比汇总
# ============================================================================

neutro_cat("【第10步】UVMR vs MVMR最终对比汇总...\n\n")

exp_labels <- c(
  neutrophil = "中性粒细胞 (ieu-b-34)",
  lymphocyte = "淋巴细胞   (ieu-b-32)",
  platelet   = "血小板     (ukb-d-30080_irnt)"
)

neutro_cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
neutro_cat(sprintf("%-30s %-6s %-6s %-14s %s\n",
            "暴露变量", "分析", "OR", "95%CI", "P值"))
neutro_cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

for (exp_name in names(harmonised_list)) {
  
  # UVMR IVW结果
  uvmr_row <- uvmr_results_list[[exp_name]]
  ivw_row  <- uvmr_row[grepl("Inverse variance weighted", uvmr_row$method), ]
  
  neutro_cat(sprintf("%-30s %-6s %-6s %-14s %s\n",
              exp_labels[exp_name], "UVMR",
              round(exp(ivw_row$b), 3),
              paste0(round(exp(ivw_row$b - 1.96*ivw_row$se), 3),
                     "–",
                     round(exp(ivw_row$b + 1.96*ivw_row$se), 3)),
              formatC(ivw_row$pval, format="e", digits=2)))
  
  # MVMR结果
  mvmr_row <- mvmr_res$result[
    grepl(exposures[[exp_name]], mvmr_res$result$id.exposure), ]
  
  if (nrow(mvmr_row) > 0) {
    neutro_cat(sprintf("%-30s %-6s %-6s %-14s %s\n",
                "", "MVMR",
                round(mvmr_row$OR, 3),
                paste0(round(mvmr_row$OR_lci, 3), "–", round(mvmr_row$OR_uci, 3)),
                formatC(mvmr_row$pval, format="e", digits=2)))
  }
  neutro_cat("─────────────────────────────────────────────────────────\n")
}

# ============================================================================
# 完成
# ============================================================================
# 提取所有MR结果

