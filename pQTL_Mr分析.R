# ============================================================
# ============================================================
rm(list = ls())

###写入Token
token <- "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiI0NzY4OTYxNDJAcXEuY29tIiwiaWF0IjoxNzg4MDgzMjg2LCJleHAiOjE3ODkyOTI4ODZ9.Fyt4awhW1fazfRlOqfqNChCzdqfnclBHHt5z8qCQ32ZSbFcoYC2R1OExZSFmbVqHVMAjgphVsmlU6HjVKUzQ9MTHoNZgPs8Bko0qjJr1MCSaudbIxBU89kv0wKsq9VstNdS7efXEGais48v4nToTHks3BUwWPdNSKKD1AbY9tcDjZ8Tp73kyA6TTb8OfrGtOYuBDfImQJhz4HpzEWiGi-fnUuBT0Z7PI9lQ2v4YFr9jaf4HHUwyynDZlRh9qlpjrplzGTFK65cngWzvCbeA87J74YT08nTX1BMGpl7PidXXzRzySAIEzPPbLf3XWUhbY2wv4uLJVEY0VOH_llefSew"              # 从网页复制的完整字符串
r_env <- file.path(Sys.getenv("HOME"), ".Renviron")
if (!file.exists(r_env)) file.create(r_env)
write(paste0("OPENGWAS_JWT=", token), r_env, append = TRUE)
message("已写入: ", r_env)
ieugwasr::get_opengwas_jwt()
remotes::install_github("mrcieu/gwasvcf")

suppressPackageStartupMessages({
  library(tidyverse)
  library(TwoSampleMR)
  library(VariantAnnotation)
  library(Rsamtools)
})
set.seed(123)
setwd("D:/Bioinfor/colorectal/MR")
mr_dir <- "D:/Bioinfor/colorectal/MR"
iv_file <- file.path(mr_dir, "IVs_selected_cis_pQTL.csv")
outcome_vcf <- file.path(mr_dir, "ebi-a-GCST90013862.vcf.gz")
deg_file <- "D:/Bioinfor/colorectal/results/GSE196006_GSE146009_GSE110224/DEG_Tumor_vs_Normal_significant.csv"
result_dir <- file.path(mr_dir, "results", "GCST90013862")
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(iv_file), file.exists(outcome_vcf),
          file.exists(paste0(outcome_vcf, ".tbi")), file.exists(deg_file))

# ---- 读取已完成的工具变量筛选结果 ----
iv_raw <- readr::read_csv(iv_file, show_col_types = FALSE) |>
  transmute(
    id.exposure, exposure, SNP,
    beta.exposure = as.numeric(beta.exposure),
    se.exposure = as.numeric(se.exposure),
    effect_allele.exposure = toupper(effect_allele.exposure),
    other_allele.exposure = toupper(other_allele.exposure),
    eaf.exposure = as.numeric(eaf.exposure),
    pval.exposure = as.numeric(pval.exposure),
    samplesize.exposure = as.numeric(samplesize.exposure),
    chr = str_remove(as.character(chr.exposure), "^chr"),
    pos = as.integer(pos.exposure)
  ) |>
  filter(!is.na(SNP), !is.na(chr), !is.na(pos),
         !is.na(beta.exposure), !is.na(se.exposure)) |>
  distinct(SNP, .keep_all = TRUE)
message("IV rows: ", nrow(iv_raw), "; proteins: ", n_distinct(iv_raw$exposure))

# ---- 使用 VCF tabix 索引提取 IV 的 outcome 关联 ----
library(ieugwasr)

gwasinfo("ebi-a-GCST90013862")  # 拿 ncase / ncontrol / samplesize / build
extract_outcome <- function(iv, vcf_file) {
  query <- GRanges(seqnames = iv$chr,
                   ranges = IRanges(start = iv$pos, end = iv$pos))
  param <- ScanVcfParam(which = query,
                        geno = c("ES", "SE", "LP", "AF"))
  vcf <- readVcf(vcf_file, genome = "hg19", param = param)
  if (!length(vcf)) return(tibble())
  tibble(
    SNP = as.character(rowRanges(vcf)$ID),
    effect_allele.outcome = toupper(as.character(unlist(alt(vcf)))),
    other_allele.outcome = toupper(as.character(ref(vcf))),
    beta.outcome = as.numeric(geno(vcf)$ES[, 1]),
    se.outcome = as.numeric(geno(vcf)$SE[, 1]),
    pval.outcome = 10^(-as.numeric(geno(vcf)$LP[, 1])),
    eaf.outcome = as.numeric(geno(vcf)$AF[, 1])
  ) |>
    filter(!is.na(SNP), !str_detect(effect_allele.outcome, ","),
           !str_detect(other_allele.outcome, ",")) |>
    distinct(SNP, .keep_all = TRUE)
}

outcome_raw <- extract_outcome_data(
  snps     = iv_raw$SNP,        # 字符向量，API 路径按 rsID 查，不受 build 影响
  outcomes = "ebi-a-GCST90013862",
  proxies  = TRUE
) |> 
  mutate(id.outcome = "ebi-a-GCST90013862",
         outcome = "Colorectal cancer",
         samplesize.outcome = 407746,
         ncase.outcome = 4161,
         ncontrol.outcome = 403585)


readr::write_csv(outcome_raw, file.path(result_dir, "CRC_outcome_for_selected_IVs.csv"))
message("Outcome SNPs extracted: ", nrow(outcome_raw))

# ---- 格式化、等位基因协调与 MR ----
names(iv_raw)

exposure_dat <- iv_raw |>
  select(id.exposure, exposure, SNP, beta.exposure, se.exposure,
         effect_allele.exposure, other_allele.exposure, eaf.exposure,
         pval.exposure, samplesize.exposure)

outcome_dat <- outcome_raw 


harmonised <- harmonise_data(exposure_dat, outcome_dat) |>
  filter(mr_keep)

if (!nrow(harmonised)) stop("IV 与 outcome 没有可协调的 SNP 交集。")
readr::write_csv(harmonised, file.path(result_dir, "harmonised_selected_IVs.csv"))

mr_results <- harmonised |>
  group_by(id.exposure, exposure, id.outcome, outcome) |>
  group_split(.keep = TRUE) |>
  purrr::map_dfr(function(dat) {
    tryCatch(
      mr(dat, method_list = c("mr_ivw", "mr_egger_regression",
                              "mr_weighted_median", "mr_weighted_mode",
                              "mr_simple_mode")),
      error = function(e) tibble()
    )
  }) |>
  left_join(harmonised |>
              count(id.exposure, exposure, id.outcome, outcome, name = "nsnp"),
            by = c("id.exposure", "exposure", "id.outcome", "outcome")) |>
  group_by(method) |>
  mutate(pval_fdr = p.adjust(pval, method = "fdr")) |>
  ungroup()
readr::write_csv(mr_results, file.path(result_dir, "MR_all_methods.csv"))

# ---- IVW 主结果与 DEG 交集 ----
ivw_results <- mr_results |>
  filter(method == "Inverse variance weighted") |>
  mutate(OR = exp(b), OR_low = exp(b - 1.96 * se),
         OR_high = exp(b + 1.96 * se),
         MR_significant = pval < 0.05)
readr::write_csv(ivw_results, file.path(result_dir, "MR_IVW_results.csv"))

deg <- readr::read_csv(deg_file, show_col_types = FALSE) |>
  transmute(GeneSymbol = toupper(str_trim(GeneSymbol)),
            DEG_logFC = logFC, DEG_FDR = adj.P.Val, DEG_Category = Category) |>
  filter(!is.na(GeneSymbol), GeneSymbol != "") |>
  distinct(GeneSymbol, .keep_all = TRUE)
mr_deg_intersection <- ivw_results |>
  filter(MR_significant) |>
  mutate(GeneSymbol = toupper(str_trim(exposure))) |>
  inner_join(deg, by = "GeneSymbol") |>
  arrange(pval_fdr)
readr::write_csv(mr_deg_intersection,
                 file.path(result_dir, "MR_IVW_FDR05_DEG_intersection.csv"))
######-------------------------韦恩图绘制--------------------------------------
library(VennDiagram)
library(dplyr)

# ---------- 1. 读取数据 ----------
mr_dir     <- "D:/Bioinfor/colorectal/MR"
deg_file   <- "D:/Bioinfor/colorectal/results/GSE196006_GSE146009_GSE110224/DEG_Tumor_vs_Normal_significant.csv"
result_dir <- file.path(mr_dir, "results", "GCST90013862")

mr_res <- read.csv(file.path(result_dir, "MR_IVW_results.csv"))
mr_sig <- mr_res %>%
  filter(MR_significant == TRUE) %>%
  pull(exposure) %>% toupper() %>% unique()

deg <- read.csv(deg_file) %>%
  mutate(GeneSymbol = toupper(trimws(GeneSymbol))) %>%
  filter(!is.na(GeneSymbol), GeneSymbol != "") %>%
  distinct(GeneSymbol, .keep_all = TRUE)

deg_up   <- deg %>% filter(Category == "Up")   %>% pull(GeneSymbol)
deg_down <- deg %>% filter(Category == "Down") %>% pull(GeneSymbol)

# ---------- 2. 公共参数（重点修改标签位置） ----------
venn_common <- list(
  scaled = FALSE,
  ext.percent = 0.70,
  ext.text = FALSE,
  
  fill   = c("#4e79a7", "#b07aa1"),
  alpha  = 0.55,
  col    = c("#4e79a7", "#b07aa1"),
  lwd    = 2.5,
  lty    = "solid",
  
  cex          = 1.3,
  fontface     = "bold",
  cat.cex      = 0.70,   # 标签进一步缩小
  cat.fontface = "bold",
  cat.col      = c("#4e79a7", "#b07aa1"),
  # =========关键参数=========
  cat.pos = c(-20, 20), # 两个分类标签角度，左圈-20°，右圈+20°，向内收，不贴画布边缘
  cat.dist = 0.08,      # 标签距离圆圈的距离，向内拉近
  margin = 0.12,        # 大幅增加四周留白！！
  print.mode = "raw",
  main.cex = 0.9,
  main.fontface = "bold"
)

# ---------- 3. 图A ----------
do.call(venn.diagram, c(
  list(
    x = list(DEG_up = deg_up, MR_significant = mr_sig),
    filename = file.path(result_dir, "Venn_A_DEGup_MRsig.png"),
    height = 12, width = 18, units = "cm", resolution = 300,
    compression = "lzw",
    main = "(A) DEG_up vs MR_significant"
  ),
  venn_common
))

# ---------- 4. 图B ----------
do.call(venn.diagram, c(
  list(
    x = list(DEG_down = deg_down, MR_significant = mr_sig),
    filename = file.path(result_dir, "Venn_B_DEGdown_MRsig.png"),
    height = 12, width = 18, units = "cm", resolution = 300,
    compression = "lzw",
    main = "(B) DEG_down vs MR_significant"
  ),
  venn_common
))

# ---------- 5. 交集基因输出 ----------
inter_up   <- intersect(deg_up, mr_sig)
inter_down <- intersect(deg_down, mr_sig)
cat("上调交集 (", length(inter_up), "):", inter_up, "\n")
cat("下调交集 (", length(inter_down), "):", inter_down, "\n")
#####-------绘制森林图-----------------------###################################
library(tidyverse)
library(patchwork)

mr_dat <- tibble(
  exposure = c("GREM1","HHIP","LAMC2","ITIH5","NTRK3","ANGPTL1","UGDH","KLB" ),
  nsnp     = c(7,2,4,6,2,2,2,4),
  method   = rep("IVW", 8),
  pval     = c(6.638436e-06, 4.769232e-05, 8.704422e-04, 2.605977e-03, 1.390666e-02, 4.471245e-02, 4.123032e-02, 3.671420e-02),
  or       = c(1.3067461, 1.7850550, 0.7587666, 0.8091456, 0.5044520, 1.1918740, 1.1953791, 1.1518175),
  or_lower = c(1.1631594, 1.3500924, 0.6449527, 0.7049400, 0.2924220, 1.0041489, 1.0071252, 1.0087617),
  or_upper = c(1.4680579, 2.3601505, 0.8926652, 0.9287552, 0.8702214,1.4146942, 1.4188218, 1.3151605)
) %>%
  mutate(
    `OR(95% CI)` = paste0(sprintf("%.3f", or), " (",
                          sprintf("%.3f", or_lower), " to ",
                          sprintf("%.3f", or_upper), ")"),
    pval_text = sprintf("%.3f", pval),
    col = "#E63946",
    y = nrow(.):1
  )

# 固定每行高度，避免森林图和文本列错位
row_height <- 0.8
y_min <- min(mr_dat$y) - 0.5
y_max <- max(mr_dat$y) + 0.5

# -------- 左侧表格 --------
p_table_left <- mr_dat %>%
  select(y, exposure, nsnp, method, pval_text) %>%
  mutate(
    nsnp = as.character(nsnp),
    pval_text = as.character(pval_text)
  ) %>%
  pivot_longer(-y, names_to = "colname", values_to = "val") %>%
  mutate(colname = factor(colname,
                          levels = c("exposure", "nsnp", "method", "pval_text"),
                          labels = c("exposure", "nsnp", "method", "pval"))) %>%
  ggplot(aes(x = colname, y = y)) +
  geom_rect(data = ~ filter(.x, y %% 2 == 0),
            aes(xmin = -Inf, xmax = Inf, ymin = y - 0.5, ymax = y + 0.5),
            fill = "#f0f0f0", inherit.aes = FALSE) +
  geom_text(aes(label = val), size = 3.2, hjust = 0.5) +
  scale_x_discrete(limits = c("exposure", "nsnp", "method", "pval")) +
  scale_y_continuous(limits = c(y_min, y_max), expand = c(0, 0)) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 9),
    axis.title = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )

# -------- 中间森林图：修复每行显示不全 --------
p_forest <- ggplot(mr_dat, aes(y = y)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = or_lower, xmax = or_upper),
                 height = 0.12, linewidth = 0.6) +
  geom_point(aes(x = or, color = col), size = 2.5) +
  scale_color_identity() +
  scale_x_continuous(trans = "log",
                     breaks = c(0.2, 0.6, 1.0, 1.4, 1.8, 2.2, 2.5),
                     limits = c(0.2, 2.6)) +
  scale_y_continuous(limits = c(y_min, y_max), expand = c(0, 0)) +
  labs(x = "Odds Ratio") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 9),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = 10),
    plot.margin = margin(0, 0, 0, 0)
  )

# -------- 右侧 OR(95% CI) 文本列 --------
p_table_right <- mr_dat %>%
  select(y, `OR(95% CI)`) %>%
  pivot_longer(-y, names_to = "colname", values_to = "val") %>%
  ggplot(aes(x = colname, y = y)) +
  geom_rect(data = ~ filter(.x, y %% 2 == 0),
            aes(xmin = -Inf, xmax = Inf, ymin = y - 0.5, ymax = y + 0.5),
            fill = "#f0f0f0", inherit.aes = FALSE) +
  geom_text(aes(label = val), size = 3.2, hjust = 0.5) +
  scale_x_discrete(limits = "OR(95% CI)") +
  scale_y_continuous(limits = c(y_min, y_max), expand = c(0, 0)) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 9),
    axis.title = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )

# -------- 拼接 --------
p_final <- wrap_plots(
  p_table_left,
  p_forest,
  p_table_right,
  ncol = 3,
  widths = c(4, 2.5, 3.5)
)

print(p_final)

ggsave("MR_forest_Fig5.png", p_final,
       width = 18, height = 7, dpi = 300)

##############################################################################

summary_tbl <- tibble(
  metric = c("IV proteins", "IV SNPs", "Outcome SNPs extracted",
             "Harmonised SNP rows", "MR proteins", "IVW FDR<0.05",
             "MR-DEG intersection"),
  value = c(n_distinct(iv_raw$exposure), n_distinct(iv_raw$SNP),
            nrow(outcome_raw), nrow(harmonised), n_distinct(mr_results$exposure),
            n_distinct(ivw_results$exposure[ivw_results$MR_significant]),
            nrow(mr_deg_intersection))
)
readr::write_csv(summary_tbl, file.path(result_dir, "MR_run_summary.csv"))
message("MR completed. Results written to: ", result_dir)
# ---- 计数 ----
cat("IVW FDR<0.05 蛋白数:", n_distinct(mr_results$exposure[mr_results$method=="Inverse variance weighted" & mr_results$pval_fdr < 0.05]), "\n")
cat("IVW raw p<0.05 蛋白数:", n_distinct(ivw_results$exposure[ivw_results$MR_significant]), "\n")
cat("MR-DEG 交集基因数:", nrow(mr_deg_intersection), "\n")

mr_deg_intersection |>
  select(exposure, pval, pval_fdr, OR, OR_low, OR_high, DEG_logFC, DEG_FDR, DEG_Category)

summary_tbl <- tibble(
  metric = c("IV proteins", "IV SNPs", "Outcome SNPs extracted",
             "Harmonised SNP rows", "MR proteins", "IVW raw p<0.05",
             "MR-DEG intersection"),
  value = c(n_distinct(iv_raw$exposure), n_distinct(iv_raw$SNP),
            nrow(outcome_raw), nrow(harmonised), n_distinct(mr_results$exposure),
            n_distinct(ivw_results$exposure[ivw_results$MR_significant]),
            nrow(mr_deg_intersection))
)
########**论文三线表**：8 基因 × IVW / MR-Egger / MR-PRESSO / 异质性 Q / 留一法稳健性
