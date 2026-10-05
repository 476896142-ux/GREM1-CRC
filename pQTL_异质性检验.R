library(dplyr)
library(ggplot2)
library(readr)
library(forestplot)

set.seed(20240918)       # MR-PRESSO 置换种子，保证可复现
options(stringsAsFactors = FALSE)
## ---------- 1. 路径配置 ----------
path_harmonised <- "results/GCST90013862/harmonised_selected_IVs.csv"
path_main_res   <- "results/GCST90013862//MR_IVW_FDR05_DEG_intersection.csv"
## ---------- 2. 读取 harmonised 数据 ----------
## harmonised CSV 已包含 TwoSampleMR 要求的全部列：
##   SNP, beta.exposure, beta.outcome, se.exposure, se.outcome,
##   eaf.exposure, exposure, id.exposure, id.outcome, outcome, ...
## 直接用 read_csv 读入后强制类型即可；如列名严格一致，也可直接喂给 mr()。
dat <- readr::read_csv(path_harmonised, show_col_types = FALSE) %>%
  as.data.frame()
## 仅保留本次 8 个目标基因的工具变量（剔除 CSV 中其他蛋白）
dat <- dat %>% filter(exposure %in% gene_order)
cat("[INFO] 载入 harmonised 数据：", nrow(dat), "行，",
    length(unique(dat$exposure)), "个暴露\n", sep = "")
## 强制 TwoSampleMR 期望的关键列为数值型
num_cols <- c("beta.exposure","beta.outcome","se.exposure","se.outcome",
              "eaf.exposure","eaf.outcome")
for (cc in num_cols) dat[[cc]] <- suppressWarnings(as.numeric(dat[[cc]]))
## ---------- 3. 主 MR：IVW + MR-Egger ----------
mr_res <- TwoSampleMR::mr(dat,
                          method_list = c("mr_ivw", "mr_egger_regression"))
View(mr_res)
cat("[INFO] mr() 完成，共", nrow(mr_res), "条结果\n")
## ---------- 4. 异质性检验（Cochran's Q）----------
het_res <- TwoSampleMR::mr_heterogeneity(dat,
                                         method_list = c("mr_ivw","mr_egger_regression"))
View(het_res)
## ---------- 5. Egger 截距（定向多效性）----------
## 仅当 nSNP >= 3 时 Egger 截距可识别；TwoSampleMR 内部会自动处理，但这里
## 显式跳过 nSNP<3 的基因以与 Python 结果保持一致。
n_snp_tab <- dat %>% group_by(exposure) %>% summarise(nSNP = n(), .groups="drop")
skip_egger <- n_snp_tab$exposure[ n_snp_tab$nSNP < 3 ]
cat("[INFO] 以下基因 SNP<3，跳过 Egger 截距：",
    paste(skip_egger, collapse = ", "), "\n")
pleio_list <- list()
for (g in gene_order) {
  if (g %in% skip_egger) {
    pleio_list[[g]] <- data.frame(
      exposure = g, id.outcome = dat$id.outcome[1],
      egger_intercept = NA, se = NA, pval = NA,
      nsnp = n_snp_tab$nSNP[match(g, n_snp_tab$exposure)]
    )
    next
  }
  sub <- dat %>% filter(exposure == g)
  pp  <- TwoSampleMR::mr_pleiotropy_test(sub)
  pleio_list[[g]] <- pp
}
pleio_res <- do.call(rbind, pleio_list)
for (g in gene_order) {
  if (g %in% skip_egger) {
    # 必须和mr_pleiotropy_test输出列完全对齐！
    pleio_list[[g]] <- data.frame(
      id.exposure = g,
      id.outcome = dat$id.outcome[1],
      outcome = dat$outcome[1],
      exposure = g,
      egger_intercept = NA,
      se = NA,
      pval = NA,
      nsnp = n_snp_tab$nSNP[match(g, n_snp_tab$exposure)] # 额外放nsnp，后面bind_rows兼容
    )
    next
  }
  sub <- dat %>% filter(exposure == g)
  pp  <- TwoSampleMR::mr_pleiotropy_test(sub)
  pp$nsnp <- n_snp_tab$nSNP[match(g, n_snp_tab$exposure)]
  pleio_list[[g]] <- pp
}
View(pleio_list)
View(pleio_list)
## ---------- 6. MR-PRESSO（仅对 nSNP >= 4 的 4 个基因）----------
## 关键参数：
##   NbDistribution = 5000   残差置换次数（R 包默认 1000，这里与 Python 一致用 5000）
##   SignifThreshold = 0.05 离群值显著性阈值
presso_list <- list()
skip_presso <- n_snp_tab$exposure[ n_snp_tab$nSNP < 4 ]
cat("[INFO] 以下基因 SNP<4，跳过 MR-PRESSO：",
    paste(skip_presso, collapse = ", "), "\n")
for (g in gene_order) {
  if (g %in% skip_presso) {
    presso_list[[g]] <- data.frame(
      exposure = g, `MR-PRESSO p-value` = NA,
      `MR-PRESSO number of outliers` = NA,
      check.names = FALSE
    )
    next
  }
  sub <- dat %>% filter(exposure == g)
  ## MRPRESSO::mr_presso() 官方签名：
  ##   mr_presso(BetaOutcome, BetaExposure, SdOutcome, SdExposure,
  ##             data = NULL, OUTLIERtest = FALSE, DISTORTIONtest = FALSE,
  ##             SignifThreshold = 0.05, NbDistribution = 1000, seed = NULL)
  ## 四个位置参数均为 data 中的列名字符串，必须先构造 data.frame 再用 data= 传入。
  presso_dat <- data.frame(
    BetaExposure = sub$beta.exposure,
    BetaOutcome  = sub$beta.outcome,
    SdExposure   = sub$se.exposure,
    SdOutcome    = sub$se.outcome      # 注意：是 outcome 的 SE
  )
  res_presso <- MRPRESSO::mr_presso(
    BetaOutcome     = "BetaOutcome",
    BetaExposure    = "BetaExposure",
    SdOutcome       = "SdOutcome",
    SdExposure      = "SdExposure",
    data            = presso_dat,
    OUTLIERtest     = TRUE,
    DISTORTIONtest  = TRUE,
    SignifThreshold = 0.05,
    NbDistribution  = 5000,   # 残差置换次数（R 包默认 1000，这里与 Python 一致用 5000）
    seed            = 20240918
  )
  ## 全局检验 p 值位于 MR-PRESSO results -> Global Test -> p.value
  global_p <- res_presso$`MR-PRESSO results`$`Global Test`$p.value
  ## 离群 SNP 数：Outlier Test 中 p < SignifThreshold 的 SNP 个数
  outlier_p <- res_presso$`MR-PRESSO results`$`Outlier Test`$p.value
  n_out     <- sum(outlier_p < 0.05, na.rm = TRUE)
  presso_list[[g]] <- data.frame(
    exposure = g,
    `MR-PRESSO p-value` = global_p,
    `MR-PRESSO number of outliers` = n_out,
    check.names = FALSE
  )
}
View(presso_dat)
set.seed(20240918)       # MR-PRESSO 置换种子，保证可复现
options(stringsAsFactors = FALSE)
## ---------- 1. 路径配置 ----------
path_harmonised <- "results/GCST90013862/harmonised_selected_IVs.csv"
path_main_res   <- "results/GCST90013862//MR_IVW_FDR05_DEG_intersection.csv"
gene_order <- c("GREM1","HHIP","LAMC2","ITIH5","NTRK3","ANGPTL1","UGDH","KLB")
## ---------- 2. 读取 harmonised 数据 ----------
## harmonised CSV 已包含 TwoSampleMR 要求的全部列：
##   SNP, beta.exposure, beta.outcome, se.exposure, se.outcome,
##   eaf.exposure, exposure, id.exposure, id.outcome, outcome, ...
## 直接用 read_csv 读入后强制类型即可；如列名严格一致，也可直接喂给 mr()。
dat <- readr::read_csv(path_harmonised, show_col_types = FALSE) %>%
  as.data.frame()
## 仅保留本次 8 个目标基因的工具变量（剔除 CSV 中其他蛋白）
dat <- dat %>% filter(exposure %in% gene_order)
cat("[INFO] 载入 harmonised 数据：", nrow(dat), "行，",
    length(unique(dat$exposure)), "个暴露\n", sep = "")
## 强制 TwoSampleMR 期望的关键列为数值型
num_cols <- c("beta.exposure","beta.outcome","se.exposure","se.outcome",
              "eaf.exposure","eaf.outcome")
for (cc in num_cols) dat[[cc]] <- suppressWarnings(as.numeric(dat[[cc]]))
## ---------- 3. 主 MR：IVW + MR-Egger ----------
mr_res <- TwoSampleMR::mr(dat,
                          method_list = c("mr_ivw", "mr_egger_regression"))
cat("[INFO] mr() 完成，共", nrow(mr_res), "条结果\n")
## ---------- 4. 异质性检验（Cochran's Q）----------
het_res <- TwoSampleMR::mr_heterogeneity(dat,
                                         method_list = c("mr_ivw","mr_egger_regression"))
## ---------- 5. Egger 截距（定向多效性）----------
## 仅当 nSNP >= 3 时 Egger 截距可识别；TwoSampleMR 内部会自动处理，但这里
## 显式跳过 nSNP<3 的基因以与 Python 结果保持一致。
n_snp_tab <- dat %>% group_by(exposure) %>% summarise(nSNP = n(), .groups="drop")
skip_egger <- n_snp_tab$exposure[ n_snp_tab$nSNP < 3 ]
cat("[INFO] 以下基因 SNP<3，跳过 Egger 截距：",
    paste(skip_egger, collapse = ", "), "\n")
pleio_list <- list()
for (g in gene_order) {
  if (g %in% skip_egger) {
    # 必须和mr_pleiotropy_test输出列完全对齐！
    pleio_list[[g]] <- data.frame(
      id.exposure = g,
      id.outcome = dat$id.outcome[1],
      outcome = dat$outcome[1],
      exposure = g,
      egger_intercept = NA,
      se = NA,
      pval = NA,
      nsnp = n_snp_tab$nSNP[match(g, n_snp_tab$exposure)] # 额外放nsnp，后面bind_rows兼容
    )
    next
  }
  sub <- dat %>% filter(exposure == g)
  pp  <- TwoSampleMR::mr_pleiotropy_test(sub)
  pp$nsnp <- n_snp_tab$nSNP[match(g, n_snp_tab$exposure)]
  pleio_list[[g]] <- pp
}
## ---------- 6. MR-PRESSO（仅对 nSNP >= 4 的 4 个基因）----------
## 关键参数：
##   NbDistribution = 5000   残差置换次数（R 包默认 1000，这里与 Python 一致用 5000）
##   SignifThreshold = 0.05 离群值显著性阈值
presso_list <- list()
skip_presso <- n_snp_tab$exposure[ n_snp_tab$nSNP < 4 ]
cat("[INFO] 以下基因 SNP<4，跳过 MR-PRESSO：",
    paste(skip_presso, collapse = ", "), "\n")
for (g in gene_order) {
  if (g %in% skip_presso) {
    presso_list[[g]] <- data.frame(
      exposure = g, `MR-PRESSO p-value` = NA,
      `MR-PRESSO number of outliers` = NA,
      check.names = FALSE
    )
    next
  }
  sub <- dat %>% filter(exposure == g)
  ## MRPRESSO::mr_presso() 官方签名：
  ##   mr_presso(BetaOutcome, BetaExposure, SdOutcome, SdExposure,
  ##             data = NULL, OUTLIERtest = FALSE, DISTORTIONtest = FALSE,
  ##             SignifThreshold = 0.05, NbDistribution = 1000, seed = NULL)
  ## 四个位置参数均为 data 中的列名字符串，必须先构造 data.frame 再用 data= 传入。
  presso_dat <- data.frame(
    BetaExposure = sub$beta.exposure,
    BetaOutcome  = sub$beta.outcome,
    SdExposure   = sub$se.exposure,
    SdOutcome    = sub$se.outcome      # 注意：是 outcome 的 SE
  )
  res_presso <- MRPRESSO::mr_presso(
    BetaOutcome     = "BetaOutcome",
    BetaExposure    = "BetaExposure",
    SdOutcome       = "SdOutcome",
    SdExposure      = "SdExposure",
    data            = presso_dat,
    OUTLIERtest     = TRUE,
    DISTORTIONtest  = TRUE,
    SignifThreshold = 0.05,
    NbDistribution  = 5000,   # 残差置换次数（R 包默认 1000，这里与 Python 一致用 5000）
    seed            = 20240918
  )
  ## 全局检验 p 值位于 MR-PRESSO results -> Global Test -> p.value
  global_p <- res_presso$`MR-PRESSO results`$`Global Test`$p.value
  ## 离群 SNP 数：Outlier Test 中 p < SignifThreshold 的 SNP 个数
  outlier_p <- res_presso$`MR-PRESSO results`$`Outlier Test`$p.value
  n_out     <- sum(outlier_p < 0.05, na.rm = TRUE)
  presso_list[[g]] <- data.frame(
    exposure = g,
    `MR-PRESSO p-value` = global_p,
    `MR-PRESSO number of outliers` = n_out,
    check.names = FALSE
  )
}
presso_res <- do.call(rbind, presso_list)
View(sub)
## ---------- 7. 留一法 ----------
loo_res <- TwoSampleMR::mr_leave_one_out(dat,
                                         method_list = c("mr_ivw"))
## ---------- 6. MR-PRESSO（仅对 nSNP >= 4 的 4 个基因）----------
## 关键参数：
##   NbDistribution = 5000   残差置换次数（R 包默认 1000，这里与 Python 一致用 5000）
##   SignifThreshold = 0.05 离群值显著性阈值
presso_list <- list()
for (g in gene_order) {
  if (g %in% skip_presso) {
    presso_list[[g]] <- data.frame(
      exposure = g,
      `MR-PRESSO p-value` = NA_real_,
      `MR-PRESSO number of outliers` = NA_integer_,
      check.names = FALSE
    )
    next
  }
  sub <- dat %>% filter(exposure == g)
  # MR-PRESSO最低要求SNP>=4，少于4个直接跳过计算
  if(nrow(sub) < 4){
    presso_list[[g]] <- data.frame(
      exposure = g,
      `MR-PRESSO p-value` = NA_real_,
      `MR-PRESSO number of outliers` = NA_integer_,
      check.names = FALSE
    )
    next
  }
  presso_dat <- data.frame(
    BetaExposure = sub$beta.exposure,
    BetaOutcome  = sub$beta.outcome,
    SdExposure   = sub$se.exposure,
    SdOutcome    = sub$se.outcome
  )
  # tryCatch捕获MR-PRESSO内部报错
  res_presso <- tryCatch({
    MRPRESSO::mr_presso(
      BetaOutcome     = "BetaOutcome",
      BetaExposure    = "BetaExposure",
      SdOutcome       = "SdOutcome",
      SdExposure      = "SdExposure",
      data            = presso_dat,
      OUTLIERtest     = TRUE,
      DISTORTIONtest  = TRUE,
      SignifThreshold = 0.05,
      NbDistribution  = 5000,
      seed            = 20240918
    )
  }, error = function(e) {
    # 计算出错返回NULL
    NULL
  })
  if(!is.null(res_presso)){
    # 安全提取，若节点不存在返回NA
    global_p <- ifelse(!is.null(res_presso$`MR-PRESSO results`$`Global Test`$p.value),
                       res_presso$`MR-PRESSO results`$`Global Test`$p.value,
                       NA_real_)
    outlier_p <- res_presso$`MR-PRESSO results`$`Outlier Test`$p.value
    n_out <- sum(outlier_p < 0.05, na.rm = TRUE)
  } else {
    global_p <- NA_real_
    n_out <- NA_integer_
  }
  presso_list[[g]] <- data.frame(
    exposure = g,
    `MR-PRESSO p-value` = global_p,
    `MR-PRESSO number of outliers` = n_out,
    check.names = FALSE
  )
}
# 合并list，推荐bind_rows
presso_res <- bind_rows(presso_list)
View(presso_res)
View(res_presso)
View(presso_res)
View(presso_res)
## ---------- 7. 留一法 ----------
loo_res <- TwoSampleMR::mr_leave_one_out(dat,
                                         method_list = c("mr_ivw"))
## ---------- 7. 留一法 ----------
loo_res <- TwoSampleMR::mr_leaveoneout(dat,
                                       method_list = c("mr_ivw"))
loo_res <- TwoSampleMR::mr_leaveoneout(dat, method = mr_ivw)
View(pp)
View(pleio_list)
View(pp)
View(pp)
View(presso_res)
## ---------- 6. MR-PRESSO（仅对 nSNP >= 4 的 4 个基因）----------
## 关键参数：
##   SignifThreshold = 0.05 离群值显著性阈值
presso_list <- list()
for (g in gene_order) {
  if (g %in% skip_presso) {
    presso_list[[g]] <- data.frame(
      exposure = g,
      `MR-PRESSO p-value` = NA_real_,
      `MR-PRESSO number of outliers` = NA_integer_,
      check.names = FALSE
    )
    next
  }
  sub <- dat %>% filter(exposure == g)
  if (nrow(sub) < 4) {                      # MR-PRESSO 至少需要 4 个 SNP
    presso_list[[g]] <- data.frame(
      exposure = g,
      `MR-PRESSO p-value` = NA_real_,
      `MR-PRESSO number of outliers` = NA_integer_,
      check.names = FALSE
    )
    next
  }
  presso_dat <- data.frame(
    BetaExposure = sub$beta.exposure,
    BetaOutcome   = sub$beta.outcome,
    SdExposure    = sub$se.exposure,
    SdOutcome     = sub$se.outcome
  )
  res_presso <- tryCatch(
    MRPRESSO::mr_presso(
      BetaOutcome     = "BetaOutcome",
      BetaExposure    = "BetaExposure",
      SdOutcome       = "SdOutcome",
      SdExposure      = "SdExposure",
      data            = presso_dat,
      OUTLIERtest     = TRUE,
      DISTORTIONtest  = TRUE,
      SignifThreshold = 0.05,
      NbDistribution  = 5000,
      seed            = 20240918
    ),
    error = function(e) NULL
  )
  if (!is.null(res_presso)) {
    # ★★★ 关键：元素名是 "Pvalue"（大写 P、大写 V），不是 "p.value"
    global_p <- res_presso[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]
    outlier_df <- res_presso[["MR-PRESSO results"]][["Outlier Test"]]
    outlier_p <- outlier_df[["Pvalue"]]     # 离群检验表的列也叫 "Pvalue"
    global_p <- if (is.null(global_p) || length(global_p) == 0) NA_real_ else global_p[1]
    n_out <- if (is.null(outlier_p) || length(outlier_p) == 0) NA_integer_
    else sum(outlier_p < 0.05, na.rm = TRUE)
  } else {
    global_p <- NA_real_
    n_out   <- NA_integer_
  }
  presso_list[[g]] <- data.frame(
    exposure = g,
    `MR-PRESSO p-value` = global_p,
    `MR-PRESSO number of outliers` = n_out,
    check.names = FALSE
  )
}
View(presso_res)
View(presso_list)
View(presso_res)
View(res_presso)
pleio_res <- do.call(rbind, pleio_list)
presso_res <- do.call(rbind, presso_list)
View(presso_res)

## ---------- 7. 留一法 ----------
loo_res <- mr_leaveoneout(dat, method = mr_ivw)
## ---------- 8. 汇总表导出 ----------
## 与 Python 版 Table 2 对应：
##   Gene, nSNP, IVW β (95%CI) P, Egger β (95%CI) P,
##   Egger 截距 estimate±SE P, Q(df) P, PRESSO 全局P/离群数, LOO 稳健性
summ <- mr_res %>%
  left_join(het_res, by = c("exposure","id.outcome","outcome","method")) %>%
  left_join(pleio_res %>% select(exposure, egger_intercept, se, pval),
            by = "exposure", suffix = c("","_int")) %>%
  left_join(presso_res, by = "exposure")
readr::write_csv(summ, file.path(out_dir, "mr_R_summary.csv"))
## ---------- 1. 路径配置 ----------
setwd("D:/Bioinfor/colorectal/MR")
mr_dir <- "D:/Bioinfor/colorectal/MR"
iv_file <- file.path(mr_dir, "IVs_selected_cis_pQTL.csv")
result_dir <- file.path(mr_dir, "results", "GCST90013862")
print(path_harmonised)
ggsave(file.path(result_dir, "leave_one_out_forest.png"),
       p_loo, width = 10, height = 8, dpi = 300)
ggsave(file.path(result_dir, "leave_one_out_forest.pdf"),
       p_loo, width = 10, height = 8)
## ---------- 8. 汇总表导出 ----------
## 与 Python 版 Table 2 对应：
##   Gene, nSNP, IVW β (95%CI) P, Egger β (95%CI) P,
##   Egger 截距 estimate±SE P, Q(df) P, PRESSO 全局P/离群数, LOO 稳健性
summ <- mr_res %>%
  left_join(het_res, by = c("exposure","id.outcome","outcome","method")) %>%
  left_join(pleio_res %>% select(exposure, egger_intercept, se, pval),
            by = "exposure", suffix = c("","_int")) %>%
  left_join(presso_res, by = "exposure")
readr::write_csv(summ, file.path(result_dir, "mr_R_summary.csv"))
cat("[INFO] 汇总表已写出：", file.path(result_dir, "mr_R_summary.csv"), "\n")
## ---------- 9. 绘图 ----------
## 9.1 漏斗图：4 个多 SNP 基因（GREM1/LAMC2/ITIH5/KLB）
for (g in c("GREM1","LAMC2","ITIH5","KLB")) {
  sub <- dat %>% filter(exposure == g)
  d_f <- sub %>%
    mutate(ratio = beta.outcome / beta.exposure,
           inv_se = 1 / se.outcome)
  ivw_b <- mr_res %>%
    filter(exposure == g, method == "MR Inverse Variance Weighted") %>%
    pull(b)
  p <- ggplot(d_f, aes(x = ratio, y = inv_se)) +
    geom_point(size = 2, color = "#2c7fb8") +
    geom_vline(xintercept = ivw_b, linetype = "dashed", color = "red") +
    labs(title = paste0("Funnel plot: ", g),
         x = "MR ratio (beta.outcome / beta.exposure)",
         y = "1 / se.outcome") +
    theme_bw(base_size = 11)
  ggsave(file.path(result_dir, paste0("funnel_", g, ".png")),
         p, width = 6, height = 5, dpi = 300)
  ggsave(file.path(result_dir, paste0("funnel_", g, ".pdf")),
         p, width = 6, height = 5)
}
## 9.2 留一法森林图（按基因分面）
## TwoSampleMR::mr_leave_one_out() 返回的剔除 SNP 列名为 removeSNP（非 rmSNP）
loo_df <- loo_res %>%
  mutate(label = ifelse(is.na(removeSNP) | removeSNP == "-",
                        "All SNPs (IVW)",
                        paste0(" -", removeSNP)))
## 9.2 留一法森林图（按基因分面）
## TwoSampleMR::mr_leave_one_out() 返回的剔除 SNP 列名为 removeSNP（非 rmSNP）
loo_df <- loo_res %>%
  mutate(label = ifelse(SNP == "All",
                        "All SNPs (IVW)",
                        paste0(" -", SNP)))
p_loo <- ggplot(loo_df, aes(x = b, y = label)) +
  geom_point(size = 1.5) +
  geom_errorbarh(aes(xmin = b - 1.96*se, xmax = b + 1.96*se), height = 0.3) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  facet_wrap(~ exposure, scales = "free_y", ncol = 2) +
  labs(title = "Leave-one-out IVW sensitivity",
       x = "β (95% CI)", y = "") +
  theme_bw(base_size = 10)
p_loo <- ggplot(loo_df, aes(x = b, y = label)) +
  geom_point(size = 1.5) +
  geom_errorbarh(aes(xmin = b - 1.96*se, xmax = b + 1.96*se), orientation = "y",height = 0.3) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  facet_wrap(~ exposure, scales = "free_y", ncol = 2) +
  labs(title = "Leave-one-out IVW sensitivity",
       x = "β (95% CI)", y = "") +
  theme_bw(base_size = 10)
ggsave(file.path(result_dir, "leave_one_out_forest.png"),
       p_loo, width = 10, height = 8, dpi = 300)
ggsave(file.path(result_dir, "leave_one_out_forest.pdf"),
       p_loo, width = 10, height = 8)
## 9.3 8 基因汇总森林图（基于 mr_res 中 IVW 固定效应）
overall <- mr_res %>%
  filter(method == "MR Inverse Variance Weighted") %>%
  mutate(lo = b - 1.96*se, hi = b + 1.96*se,
         p_text = ifelse(pval < 0.001, "<0.001", sprintf("%.3f", pval))) %>%
  arrange(match(exposure, gene_order))
p_ov <- ggplot(overall, aes(x = b, y = exposure)) +
  geom_point(size = 3, color = "#d95f0e") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.3) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  labs(title = "IVW fixed-effect estimates across 8 blood protein pQTLs",
       x = "β (95% CI)", y = "") +
  theme_bw(base_size = 11)
ggsave(file.path(result_dir, "mr_overall_forest.png"),
       p_ov, width = 7, height = 4.5, dpi = 300)
ggsave(file.path(result_dir, "mr_overall_forest.pdf"),
       p_ov, width = 7, height = 4.5)
cat("[DONE] 全部结果与图形已输出至：", result_dir, "\n")
cat("[DONE] 全部结果与图形已输出至：", result_dir, "\n")
