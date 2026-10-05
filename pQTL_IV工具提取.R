rm(list=ls())
library(tidyverse)
library(ieugwasr)
setwd("D:/Bioinfor/colorectal/MR/")
# ========== 参数配置 ==========
ld_ref <- "D:/Bioinfor/colorectal/MR/g1000_eur/g1000_eur"  # 面板前缀
plink_path <- "D:/Bioinfor/colorectal/MR/LD/plink_win64/plink.exe"

p_thresh <- 5e-8      # 顺式pQTL显著性阈值
r2_thresh <- 0.001    # LD r²阈值
kb_window <- 10000    # 10Mb窗口
f_thresh <- 10        # 弱工具F统计量阈值

out_rds  <- "pQTL_exp_dat_clump.rds"   # RDS 保存路径(保留所有列)
out_csv  <- "pQTL_exp_dat_clump.csv"   # CSV 保存路径(便于查看/共享)
# ========== 1. 读取并格式化数据 ==========
exp_dat <- read_csv("cis.exposure_data.csv", show_col_types = FALSE) %>%
  mutate(
    pval = as.numeric(pval.exposure),
    beta = as.numeric(beta.exposure),
    se = as.numeric(se.exposure),
    F_stat = (beta / se)^2
  )

# ========== 2. 按蛋白分组做独立LD clump ==========
clump_one_protein <- function(df) {
  tryCatch({
    # 关键：重命名为ieugwasr要求的标准列名 rsid + pval
    snp_input <- df %>% select(rsid = SNP, pval)
    
    clump_res <- ieugwasr::ld_clump(
      dat = snp_input,
      clump_kb = kb_window,
      clump_r2 = r2_thresh,
      clump_p = p_thresh,
      bfile = ld_ref,
      plink_bin = plink_path
    )
    
    # 返回保留的独立SNP对应行
    df %>% filter(SNP %in% clump_res$rsid)
  }, error = function(e) {
    message("蛋白 ", unique(df$exposure), " 处理失败: ", e$message)
    return(tibble())
  })
}

# 先按P值初筛，再按蛋白分组clump
exp_dat_clump <- exp_dat %>%
  filter(pval < p_thresh) %>%       # 先初筛显著SNP，减少计算量
  group_by(id.exposure, exposure) %>%
  group_split() %>%
  lapply(clump_one_protein) %>%
  bind_rows()
# ---------- 保存 ----------
saveRDS(exp_dat_clump, file = out_rds)
write.csv(exp_dat_clump, file = out_csv, row.names = FALSE)

# ========== 3. 弱工具筛选（F > 10） ==========
exp_dat_final <- exp_dat_clump %>%
  filter(F_stat > f_thresh)

# ========== 4. 结果汇总与导出 ==========
iv_summary <- exp_dat_final %>%
  group_by(id.exposure, exposure) %>%
  summarise(
    n_iv = n(),
    min_pval = min(pval),
    max_F = max(F_stat),
    chr = unique(chr.exposure),
    .groups = "drop"
  )

write_csv(exp_dat_final, "IVs_selected_cis_pQTL.csv")
write_csv(iv_summary, "IV_summary_per_protein.csv")

cat("筛选完成：共保留", nrow(iv_summary), "个蛋白，", nrow(exp_dat_final), "个独立工具变量\n")
