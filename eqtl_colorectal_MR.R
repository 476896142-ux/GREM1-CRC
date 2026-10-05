setwd("D:/Bioinfor/colorectal/MR")


#MR分析 IVW为主，MR‑Egger,weighted median,simple mode,weighted mode
mr_res <- mr(harmonised_data, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median","mr_simple_mode","mr_weighted_mode"))
head(mr_res)
write.table(mr_res, file.path("D:/Bioinfor/colorectal/MR", "mr_res.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

save(mr_res,file = "mr_res.RData")
load("mr_res.RData")
names(mr_res)
# ===== 加载依赖 =====
library(dplyr)
library(readr)
library(AnnotationDbi)
library(org.Hs.eg.db)

# ===== 1) 读取 DEG =====
deg_file <- "D:/Bioinfor/colorectal/results/GSE196006_GSE146009_GSE110224/DEG_Tumor_vs_Normal_significant.csv"
deg <- readr::read_csv(deg_file, show_col_types = FALSE) |>
  transmute(GeneSymbol = toupper(str_trim(GeneSymbol)),
            DEG_logFC = logFC, DEG_FDR = adj.P.Val, DEG_Category = Category) |>
  filter(!is.na(GeneSymbol), GeneSymbol != "") |>
  distinct(GeneSymbol, .keep_all = TRUE)

# ===== 2) 构建 gene_map =====
# 先给 mr_res 加一列去版本号的暴露 ID，后续 join 统一用它
mr_res$exposure_clean <- sub("\\..*$", "", mr_res$exposure)

gene_ids <- unique(mr_res$exposure_clean)
gene_map <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys    = gene_ids,
                                  keytype = "ENSEMBL",
                                  columns = c("SYMBOL", "GENENAME")) |>
  filter(!is.na(SYMBOL) & !is.na(GENENAME)) |> 
  distinct(ENSEMBL, .keep_all = TRUE)

# ===== 3) MR 结果：不做 FDR，直接用 pval < 0.05 筛选 =====
mr_sig <- mr_res %>%
  filter(pval < 0.05)

# ===== 4) 注释：用 exposure_clean 与 gene_map 匹配 =====
mr_sig_annot <- mr_sig %>%
  left_join(gene_map, by = c("exposure_clean" = "ENSEMBL")) %>%
  mutate(GeneSymbol = toupper(coalesce(SYMBOL, exposure_clean))) %>%
  filter(!is.na(GeneSymbol), GeneSymbol != "") %>%
  distinct(GeneSymbol, method, exposure_clean, .keep_all = TRUE)

head(mr_sig_annot)   # 此时 SYMBOL / GENENAME / GeneSymbol 应显示基因符号
# ===== 5) 与 DEG 取交集（保留双边信息） =====
names(mr_sig_annot)
ivw_eqtl_results <- mr_sig_annot |>
  filter(method == "Inverse variance weighted") |>
  mutate(OR = exp(b), OR_low = exp(b - 1.96 * se),
         OR_high = exp(b + 1.96 * se),
         MR_significant = pval < 0.05)

eqtl_deg_intersection <- ivw_eqtl_results |>
  filter(MR_significant) |>
  inner_join(deg, by = "GeneSymbol")

readr::write_csv(eqtl_deg_intersection, file.path("results", "eqtl_deg_intersection.csv"))
##########---------韦恩图-------------------------#############################
library(VennDiagram)
library(dplyr)
# ---------- 1. 读取数据 ----------
deg_file   <- "D:/Bioinfor/colorectal/results/GSE196006_GSE146009_GSE110224/DEG_Tumor_vs_Normal_significant.csv"
result_dir <- "D:/Bioinfor/colorectal/results"

mr_res <- read.csv(file.path("results", "eqtl_deg_intersection.csv"))
mr_sig <- mr_res %>%
  filter(MR_significant == TRUE) %>%
  pull(GeneSymbol) %>% toupper() %>% unique()

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
print(mr_res$OR_high)
mr_dat <- tibble(
  exposure = c("CSE1L","ZNF239","KDELR3","LAMC2","FPR1","EPHB1","WDR72","EPB41L3","MYO15B","GPER1","FABP2","MS4A8"),
  nsnp     = c(2,4,2,2,3,3,2,4,2,3,2,2),
  method   = rep("IVW", 12),
  pval     = c(0.0023951966, 0.0431908240, 0.0248061900, 0.0478810707, 0.0187892830, 0.0082401633, 0.0192900218, 0.0339193943,0.0408254440, 0.0315525149, 0.0003294808, 0.0126803729),
  or       = c(0.7820784, 1.0524239, 1.0677817, 0.9162047, 0.8844206, 0.9463075, 0.7090681, 0.8244231, 1.1566779, 1.0761603,0.9039688, 1.1639122),
  or_lower = c(0.6673296, 1.0015645, 1.0083452, 0.8401155, 0.7982963, 0.9083475, 0.5316388, 0.6897095, 1.0060887, 1.0065086,0.8555025, 1.0329655),
  or_upper = c(0.9165584, 1.1058659, 1.1307217, 0.9991854, 0.9798365, 0.9858539, 0.9457128, 0.9854489, 1.3298070, 1.1506320,0.9551810, 1.3114588)
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
  geom_errorbar(aes(xmin = or_lower, xmax = or_upper),
                height = 0.12, linewidth = 0.6, orientation = "y") +
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

















################################################################################


# ===== 0) 12 个交集基因 =====
sig_genes <- unique(eqtl_deg_intersection$GeneSymbol)
length(sig_genes)   # 应为 12

# ===== 1) 从 harmonised_data 提取这 12 个基因的数据 =====
# id.exposure 形如 "Colorectal eQTL: ENSG00000168491.10"，先解析出 ENSEMBL 再映射
harm <- harmonised_data %>%
  mutate(exposure_clean = sub("\\..*$", "",
                              sub(".*eQTL:\\s*", "", id.exposure))) %>%
  left_join(gene_map, by = c("exposure_clean" = "ENSEMBL")) %>%
  mutate(GeneSymbol = toupper(SYMBOL)) %>%
  filter(GeneSymbol %in% sig_genes)

# ===== 2) 敏感性分析 =====
# 2.1 异质性检验（Cochran's Q）：Q_pval > 0.05 → SNP 间效应一致，结果稳健
het <- mr_heterogeneity(harm)

# 2.2 水平多效性检验（MR-Egger 截距）：pval > 0.05 → 不存在水平多效性
pleio <- mr_pleiotropy_test(harm)

# 2.3 留一法：逐个剔除 SNP，看因果估计是否翻转/消失
loo <- mr_leaveoneout(harm)
mr_leaveoneout_plot(loo)   # 森林图，点不跨 0 且波动小 = 稳健

# ===== 3) Steiger 方向性检验 =====
# correct_causal_direction = TRUE → 支持"基因表达 → 结直肠癌"方向
steiger <- directionality_test(harm)

# ===== 4) 合并成一张验证汇总表 =====
id2gene <- harm %>% distinct(exposure, GeneSymbol)

map_gene <- function(d) {
  d %>% left_join(id2gene, by = "exposure") %>% select(-exposure)
}

validation_summary <- map_gene(het) %>%
  filter(method == "Inverse variance weighted") %>%
  select(GeneSymbol, Q, Q_pval) %>%
  left_join(map_gene(pleio) %>% select(GeneSymbol, Egger_intercept = egger_intercept, Egger_pval = pval),
            by = "GeneSymbol") %>%
  left_join(map_gene(steiger) %>% select(GeneSymbol, correct_direction = correct_causal_direction, Steiger_pval = steiger_pval),
            by = "GeneSymbol")
head(validation_summary)


# 保存

write.csv(validation_summary, "D:/Bioinfor/colorectal/MR/results/validation_summary.csv", row.names = FALSE)
