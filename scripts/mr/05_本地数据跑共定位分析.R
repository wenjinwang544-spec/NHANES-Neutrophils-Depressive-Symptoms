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
neutro_domain_entry('scripts/mr/05_本地数据跑共定位分析.R')
# Derived candidate; see provenance/script_selection.csv. Invoke via run_job.py.
# ============================================================================
# 共定位分析脚本 v5（seqminer版，无需Bioconductor和tabix）
# 暴露：中性粒细胞计数 ieu-b-34.vcf.gz（本地文件）
# 结局：PGC MDD2025 noUKBB（本地文件）
# ============================================================================


# ── 第0步：安装缺少的包（只需一次）──────────────────────────────────────────
# install.packages("seqminer", repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
# install.packages("coloc",    repos="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")


# ── 第1步：加载包 ─────────────────────────────────────────────────────────────
library(seqminer)
library(coloc)
library(data.table)
library(dplyr)

neutro_cat("✓ 包加载完成\n\n")


# ── 第2步：设置文件路径 ───────────────────────────────────────────────────────
NEU_VCF         <- "data/gwas/ieu-b-34.vcf.gz"
MDD_FILE        <- "data/gwas/pgc-mdd2025_no23andMe-noUKBB_eur_v3-49-24-11_formatted.tsv.gz"
HARMONISED_FILE <- "data/gwas_cache/harmonised_list.rds"
OUT_DIR         <- "work/coloc/"; neutro_dir_create(OUT_DIR,recursive=TRUE,showWarnings=FALSE)

for (f in c(NEU_VCF, MDD_FILE, HARMONISED_FILE)) {
  if (!neutro_file_exists(f)) stop(paste("找不到文件:", f))
}
neutro_cat("✓ 文件路径验证通过\n\n")


# ── 第3步：定义VCF解析函数 ────────────────────────────────────────────────────
# seqminer读出的最后一列格式：ES:SE:LP:AF:SS:ID
# ES=beta, SE=se, LP=-log10(p), AF=eaf, SS=样本量

query_neu_vcf <- function(vcf_path, chr, start, end) {

  region <- paste0(chr, ":", start, "-", end)

  raw <- tryCatch(
    neutro_tabix_table(tabixFile  = vcf_path,
                     tabixRange = region),
    error = function(e) { message("读取失败: ", e$message); NULL }
  )

  if (is.null(raw) || nrow(raw) == 0) return(NULL)

  # 解析最后一列（包含实际数值）
  val_col <- raw[[ncol(raw)]]
  vals    <- strsplit(as.character(val_col), ":")

  # 提取每个字段
  beta_vec <- as.numeric(sapply(vals, function(x) x[1]))  # ES
  se_vec   <- as.numeric(sapply(vals, function(x) x[2]))  # SE
  lp_vec   <- as.numeric(sapply(vals, function(x) x[3]))  # LP (-log10 p)
  af_vec   <- as.numeric(sapply(vals, function(x) x[4]))  # AF
  ss_vec   <- as.numeric(sapply(vals, function(x) x[5]))  # SS

  result <- data.frame(
    rsid = as.character(raw$ID),
    chr  = as.character(raw$CHROM),
    pos  = as.integer(raw$POS),
    ea   = as.character(raw$ALT),
    nea  = as.character(raw$REF),
    beta = beta_vec,
    se   = se_vec,
    pval = 10^(-lp_vec),
    eaf  = af_vec,
    n    = ss_vec,
    stringsAsFactors = FALSE
  )

  # 过滤无效行
  result <- result[!is.na(result$beta) &
                   !is.na(result$se)   &
                   result$se > 0       &
                   !is.na(result$eaf), ]
  return(result)
}

# 测试函数
neutro_cat("测试VCF读取...\n")
test <- query_neu_vcf(NEU_VCF, "1", 1000000, 1100000)
if (is.null(test) || nrow(test) == 0) stop("VCF测试失败，请检查文件")
neutro_cat("✓ VCF可读，测试区域SNP数:", nrow(test), "\n\n")


# ── 第4步：加载MDD数据 ────────────────────────────────────────────────────────
neutro_cat("读取MDD数据（1-2分钟）...\n")
mdd_full <- neutro_fread(MDD_FILE, showProgress = TRUE)
mdd_full <- mdd_full[standard_error > 0 &
                     !is.na(beta)          &
                     !is.na(standard_error)]
neutro_cat("✓ MDD有效SNP数:", nrow(mdd_full), "\n\n")


# ── 第5步：加载MR结果，提取前50个lead SNPs ───────────────────────────────────
neutro_cat("加载MR数据...\n")
harmonised_list <- neutro_readRDS(HARMONISED_FILE)
neu_harm        <- harmonised_list[["neutrophil"]]

lead_snps <- neu_harm %>%
  filter(pval.exposure < 5e-8) %>%
  arrange(pval.exposure) %>%
  slice_head(n = 50) %>%
  select(SNP, beta.exposure, se.exposure,
         beta.outcome, se.outcome, pval.exposure) %>%
  distinct()

neutro_cat("✓ 使用前", nrow(lead_snps), "个最显著lead SNP（共有",
    sum(neu_harm$pval.exposure < 5e-8), "个）\n\n")


# ── 第6步：获取lead SNP位置并合并重叠窗口 ────────────────────────────────────
neutro_cat("获取SNP位置信息...\n")

WINDOW_BP <- 500000

snp_positions <- mdd_full %>%
  filter(rsid %in% lead_snps$SNP) %>%
  select(SNP = rsid, chr = chromosome, pos = base_pair_location) %>%
  distinct() %>%
  inner_join(lead_snps %>% select(SNP, pval.exposure), by = "SNP") %>%
  arrange(chr, pos, pval.exposure)

missing_snps <- setdiff(lead_snps$SNP, snp_positions$SNP)
if (length(missing_snps) > 0) {
  neutro_cat("⚠ 以下SNP在MDD文件中未找到位置，将跳过：",
      paste(missing_snps, collapse=", "), "\n")
}

collapse_overlapping_windows <- function(x, flank_bp) {
  x <- x %>% arrange(chr, pos, pval.exposure)
  region_rows <- list()
  region_id <- 0L

  for (chr_i in unique(x$chr)) {
    z <- x[x$chr == chr_i, , drop = FALSE]
    current_members <- 1L
    current_start <- max(1, z$pos[1] - flank_bp)
    current_end <- z$pos[1] + flank_bp

    finalize_region <- function(member_idx, start_bp, end_bp) {
      members <- z[member_idx, , drop = FALSE]
      representative <- members[which.min(members$pval.exposure), , drop = FALSE]
      data.frame(
        SNP = representative$SNP,
        chr = representative$chr,
        pos = representative$pos,
        region_start = start_bp,
        region_end = end_bp,
        N_leads = nrow(members),
        Lead_variants = paste(members$SNP, collapse = ";"),
        stringsAsFactors = FALSE
      )
    }

    if (nrow(z) > 1) {
      for (j in 2:nrow(z)) {
        next_start <- max(1, z$pos[j] - flank_bp)
        next_end <- z$pos[j] + flank_bp
        if (next_start <= current_end) {
          current_members <- c(current_members, j)
          current_end <- max(current_end, next_end)
        } else {
          region_id <- region_id + 1L
          region_rows[[region_id]] <- finalize_region(
            current_members, current_start, current_end
          )
          current_members <- j
          current_start <- next_start
          current_end <- next_end
        }
      }
    }

    region_id <- region_id + 1L
    region_rows[[region_id]] <- finalize_region(
      current_members, current_start, current_end
    )
  }

  bind_rows(region_rows) %>% arrange(chr, region_start)
}

analysis_regions <- collapse_overlapping_windows(snp_positions, WINDOW_BP)

neutro_cat("✓ 获取到", nrow(snp_positions), "/", nrow(lead_snps), "个SNP位置\n")
neutro_cat("✓ 重叠窗口合并后共", nrow(analysis_regions), "个唯一基因组区域\n\n")
print(analysis_regions)
neutro_cat("\n")


# ── 第7步：共定位分析主循环 ───────────────────────────────────────────────────
neutro_cat("开始共定位分析（共", nrow(analysis_regions), "个区域）...\n\n")

MIN_SNPS  <- 50
COLOC_P1  <- 1e-4
COLOC_P2  <- 1e-4
COLOC_P12 <- 1e-5
NEU_N     <- 563946

MDD_NCASES    <- median(mdd_full$ncases,    na.rm = TRUE)
MDD_NCONTROLS <- median(mdd_full$ncontrols, na.rm = TRUE)
MDD_S         <- MDD_NCASES / (MDD_NCASES + MDD_NCONTROLS)

neutro_cat(sprintf("MDD样本：%.0f例，%.0f对照，s=%.4f\n\n",
            MDD_NCASES, MDD_NCONTROLS, MDD_S))

coloc_results <- list()

for (i in seq_len(nrow(analysis_regions))) {

  snp_id <- analysis_regions$SNP[i]
  chr    <- analysis_regions$chr[i]
  pos    <- as.numeric(analysis_regions$pos[i])

  region_start <- analysis_regions$region_start[i]
  region_end   <- analysis_regions$region_end[i]

  neutro_cat(sprintf("── 区域 %d/%d：%s  chr%s:%d-%d ──\n",
              i, nrow(analysis_regions), snp_id, chr,
              region_start, region_end))

  # 7a：提取MDD区域数据
  mdd_region <- mdd_full %>%
    filter(chromosome        == chr,
           base_pair_location >= region_start,
           base_pair_location <= region_end)

  neutro_cat("  MDD区域SNP数:", nrow(mdd_region), "\n")

  if (nrow(mdd_region) < MIN_SNPS) {
    neutro_cat("  ⚠ MDD SNP不足", MIN_SNPS, "，跳过\n\n")
    next
  }

  # 7b：从本地VCF读取中性粒细胞数据
  tryCatch({

    neu_region <- query_neu_vcf(NEU_VCF, chr, region_start, region_end)

    if (is.null(neu_region) || nrow(neu_region) == 0) {
      neutro_cat("  ⚠ 该区域无中性粒细胞数据，跳过\n\n")
      next
    }

    neutro_cat("  中性粒细胞区域SNP数:", nrow(neu_region), "\n")

    # 7c：对齐数据集
    common <- intersect(mdd_region$rsid, neu_region$rsid)
    neutro_cat("  共同SNP数:", length(common), "\n")

    if (length(common) < MIN_SNPS) {
      neutro_cat("  ⚠ 共同SNP不足", MIN_SNPS, "，跳过\n\n")
      next
    }

    mdd_sub <- mdd_region[mdd_region$rsid %in% common, ]
    neu_sub <- neu_region[neu_region$rsid %in% common, ]

    mdd_sub <- mdd_sub[order(mdd_sub$rsid), ]
    neu_sub <- neu_sub[order(neu_sub$rsid), ]

    ok <- mdd_sub$standard_error > 0 & !is.na(mdd_sub$standard_error)
    mdd_sub <- mdd_sub[ok, ]
    neu_sub <- neu_sub[ok, ]

    neutro_cat("  有效SNP数:", nrow(neu_sub), "\n")

    if (nrow(neu_sub) < MIN_SNPS) {
      neutro_cat("  ⚠ 过滤后不足", MIN_SNPS, "，跳过\n\n")
      next
    }

    # 7d：运行coloc.abf()
    result <- coloc.abf(
      dataset1 = list(
        beta    = neu_sub$beta,
        varbeta = neu_sub$se^2,
        N       = NEU_N,
        MAF     = pmin(neu_sub$eaf, 1 - neu_sub$eaf),
        type    = "quant",
        snp     = neu_sub$rsid
      ),
      dataset2 = list(
        beta    = mdd_sub$beta,
        varbeta = mdd_sub$standard_error^2,
        N       = MDD_NCASES + MDD_NCONTROLS,
        type    = "cc",
        s       = MDD_S,
        snp     = mdd_sub$rsid
      ),
      p1 = COLOC_P1, p2 = COLOC_P2, p12 = COLOC_P12
    )

    pp <- result$summary

    neutro_cat(sprintf("  PP.H0=%.3f | H1=%.3f | H2=%.3f | H3=%.3f | H4=%.3f\n",
                pp["PP.H0.abf"], pp["PP.H1.abf"], pp["PP.H2.abf"],
                pp["PP.H3.abf"], pp["PP.H4.abf"]))

    if      (pp["PP.H4.abf"] >= 0.8) neutro_cat("  ✓✓ 强共定位 PP.H4 =", round(pp["PP.H4.abf"],3), "\n")
    else if (pp["PP.H4.abf"] >= 0.5) neutro_cat("  ✓  共定位阳性 PP.H4 =", round(pp["PP.H4.abf"],3), "\n")
    else if (pp["PP.H3.abf"] >= 0.5) neutro_cat("  →  两性状独立信号 PP.H3 =", round(pp["PP.H3.abf"],3), "\n")
    else                              neutro_cat("  →  证据不足\n")

    coloc_results[[snp_id]] <- data.frame(
      Lead_SNP  = snp_id,
      CHR       = chr,
      POS       = pos,
      Region_Start = region_start,
      Region_End   = region_end,
      N_Lead_Variants = analysis_regions$N_leads[i],
      Lead_Variants = analysis_regions$Lead_variants[i],
      N_SNPs    = nrow(neu_sub),
      PP.H0     = round(pp["PP.H0.abf"], 4),
      PP.H1     = round(pp["PP.H1.abf"], 4),
      PP.H2     = round(pp["PP.H2.abf"], 4),
      PP.H3     = round(pp["PP.H3.abf"], 4),
      PP.H4     = round(pp["PP.H4.abf"], 4),
      Coloc_pos = pp["PP.H4.abf"] >= 0.5,
      Coloc_str = pp["PP.H4.abf"] >= 0.8,
      row.names = NULL
    )

  }, error = function(e) neutro_cat("  ✗ 错误:", e$message, "\n"))

  neutro_cat("\n")
}


# ── 第8步：汇总结果 ───────────────────────────────────────────────────────────
neutro_cat("══════════════════════════════════════════\n")
neutro_cat("共定位分析完成\n")
neutro_cat("══════════════════════════════════════════\n\n")

if (length(coloc_results) == 0) {

  neutro_cat("⚠ 没有成功完成的区域\n")

} else {

  summary_df           <- do.call(rbind, coloc_results)
  rownames(summary_df) <- NULL

  neutro_cat("结果汇总：\n")
  print(summary_df)
  neutro_cat("\n")
  neutro_cat("分析区域数:           ", nrow(summary_df), "\n")
  neutro_cat("PP.H4 ≥ 0.5（阳性）: ", sum(summary_df$Coloc_pos, na.rm=TRUE), "\n")
  neutro_cat("PP.H4 ≥ 0.8（强证）: ", sum(summary_df$Coloc_str, na.rm=TRUE), "\n")

  unique_region_output <- file.path(
    OUT_DIR, "coloc_results_unique_regions_N563946.csv"
  )
  neutro_write_csv(summary_df,
            unique_region_output,
            row.names = FALSE)
  neutro_cat("\n✓ 结果已保存至:", unique_region_output, "\n")

  # 文章报告参考语句
  n_pos   <- sum(summary_df$Coloc_pos, na.rm=TRUE)
  n_total <- nrow(summary_df)

  neutro_cat("\n── 文章报告参考语句 ──\n")
  if (n_pos > 0) {
    best <- summary_df[which.max(summary_df$PP.H4), ]
    neutro_cat(sprintf(
      "Colocalization analysis across the %d most significant neutrophil\n", n_total))
    neutro_cat(sprintf(
      "count loci identified %d locus/loci with evidence of shared genetic\n", n_pos))
    neutro_cat(sprintf(
      "architecture with MDD (PP.H4 ≥ 0.5), including %s (PP.H4 = %.3f),\n",
      best$Lead_SNP, best$PP.H4))
    neutro_cat("supporting shared causal variants underlying both neutrophil\n")
    neutro_cat("activation and depression risk.\n")
  } else {
    neutro_cat(sprintf(
      "Colocalization analysis across the %d most significant neutrophil\n", n_total))
    neutro_cat("count loci found no locus with strong evidence of shared causal\n")
    neutro_cat("variants with MDD (all PP.H4 < 0.5), consistent with polygenic\n")
    neutro_cat("effects operating through distinct variants and with the attenuated\n")
    neutro_cat("IVW estimate in the primary MR analysis.\n")
  }
}

neutro_cat("\n✓ 脚本运行完毕\n")
