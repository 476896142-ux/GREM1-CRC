#============================================================
# eQTL‑MR分析：12个MR显著基因，结直肠癌结局 ebi‑a‑GCST90013862
#============================================================
rm(list = ls())
library(dplyr)
library(purrr)
library(ggplot2)

out_dir <- "D:/Bioinfor/colorectal/MR/eQTL"
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

#==================== 1.提取12个MR显著基因的harmonised数据集 ====================
load("D:/Bioinfor/colorectal/MR/eQTL/eqtl_result.RData")
mr_res <- read.csv("D:/Bioinfor/colorectal/MR/results/eqtl_deg_intersection.csv")

mr_sig_symbols <- mr_res %>%
  filter(MR_significant == TRUE) %>%
  pull(GeneSymbol) %>%
  toupper() %>%
  unique()
cat("12个MR显著基因:\n")
print(mr_sig_symbols)

# ENSG ID去除版本号
ensg2symbol <- mr_res %>%
  dplyr::select(exposure, GeneSymbol) %>%
  distinct() %>%
  mutate(ensg_clean = sub("\\.[0-9]+$", "", exposure))

# 从id.exposure字符串提取纯ENSG编号
harmonised_data <- harmonised_data %>%
  mutate(
    raw_ensg = sub("^.*:\\s*", "", id.exposure),
    ensg_clean = sub("\\.[0-9]+$", "", raw_ensg)
  )

harmonised_anno <- harmonised_data %>%
  left_join(ensg2symbol, by = "ensg_clean")

harmonised_12gene <- harmonised_anno %>%
  filter(GeneSymbol %in% mr_sig_symbols, mr_keep == TRUE)

table_snp_count <- harmonised_12gene %>%
  count(GeneSymbol, name = "n_snp") %>%
  arrange(desc(n_snp))
cat("\n====每个基因可用SNP数量(mr_keep=TRUE)====\n")
print(table_snp_count)

save(harmonised_12gene, table_snp_count,
     file = file.path(out_dir, "harmonised_12gene.RData"))
write.csv(harmonised_12gene,
          file = file.path(out_dir, "harmonised_12gene.csv"), row.names = FALSE)

#==================== 2.单基因循环MR、异质性、多效性、留一法 ====================
library(TwoSampleMR)
library(MRPRESSO)

load(file.path(out_dir, "harmonised_12gene.RData"))

harmonised_clean <- harmonised_12gene %>%
  filter(
    id.outcome == "ebi-a-GCST90013862",
    !is.na(id.exposure),
    !is.na(id.outcome),
    mr_keep == TRUE
  )

harmonised_fixed <- harmonised_clean %>%
  mutate(
    GeneSymbol = toupper(GeneSymbol),
    id.exposure = GeneSymbol,
    exposure = GeneSymbol,
    id.outcome = "ebi-a-GCST90013862",
    outcome = "Colorectal cancer"
  ) %>%
  distinct(GeneSymbol, SNP, .keep_all = TRUE)

harmonised_fixed %>%
  count(GeneSymbol, id.exposure, id.outcome, outcome) %>%
  print()

# 单基因MR封装函数，带异常捕获
run_gene_mr <- function(dat_g) {
  tryCatch({
    gene_name <- unique(dat_g$GeneSymbol)
    n_snp <- nrow(dat_g)
    
    stopifnot(
      length(gene_name) == 1L,
      length(unique(dat_g$id.exposure)) == 1L,
      length(unique(dat_g$id.outcome)) == 1L,
      all(c("exposure", "outcome") %in% names(dat_g))
    )
    
    method_list <- if (n_snp == 1L) {
      "mr_wald_ratio"
    } else if (n_snp == 2L) {
      "mr_ivw"
    } else {
      c("mr_ivw", "mr_egger_regression", "mr_weighted_median")
    }
    
    mr_result <- mr(dat_g, method_list = method_list) %>%
      mutate(GeneSymbol = gene_name, n_snp_available = n_snp)
    
    # --------异质性检验 Cochran's Q--------
    heterogeneity_result <- if (n_snp >= 2L) {
      het_methods <- if (n_snp >= 3L) c("mr_ivw", "mr_egger_regression") else "mr_ivw"
      mr_heterogeneity(dat_g, method_list = het_methods) %>%
        mutate(GeneSymbol = gene_name, n_snp_available = n_snp)
    } else {
      tibble()
    }
    
    # --------MR‑Egger截距（水平多效性）--------
    pleiotropy_result <- if (n_snp >= 3L) {
      mr_pleiotropy_test(dat_g) %>%
        mutate(GeneSymbol = gene_name, n_snp_available = n_snp)
    } else {
      tibble()
    }
    
    # --------留一法敏感性分析--------
    leave_one_out_result <- if (n_snp >= 3L) {
      mr_leaveoneout(dat_g, method = mr_ivw) %>%
        mutate(GeneSymbol = gene_name, n_snp_available = n_snp)
    } else {
      tibble()
    }
    
    return(list(
      mr = mr_result,
      heterogeneity = heterogeneity_result,
      pleiotropy = pleiotropy_result,
      leave_one_out = leave_one_out_result,
      error = NULL
    ))
  }, error = function(e) {
    message(sprintf("!!! Gene %s failed: %s", unique(dat_g$GeneSymbol), e$message))
    return(list(mr = tibble(), heterogeneity = tibble(), pleiotropy = tibble(),
                leave_one_out = tibble(), error = e$message))
  })
}

analysis_by_gene <- harmonised_fixed %>%
  split(.$GeneSymbol) %>%
  map(run_gene_mr)

mr_res_12    <- map_dfr(analysis_by_gene, "mr")
mr_het_12    <- map_dfr(analysis_by_gene, "heterogeneity")
mr_pleio_12  <- map_dfr(analysis_by_gene, "pleiotropy")
mr_loo_12    <- map_dfr(analysis_by_gene, "leave_one_out")

# 输出MR全套表格
write.csv(mr_res_12,    file.path(out_dir, "mr_res_12gene.csv"), row.names = FALSE)
write.csv(mr_het_12,    file.path(out_dir, "mr_heterogeneity_12gene.csv"), row.names = FALSE)
write.csv(mr_pleio_12,  file.path(out_dir, "mr_pleiotropy_12gene.csv"), row.names = FALSE)
write.csv(mr_loo_12,    file.path(out_dir, "mr_leaveoneout_12gene.csv"), row.names = FALSE)

save(analysis_by_gene, mr_res_12, mr_het_12, mr_pleio_12, mr_loo_12,
     file = file.path(out_dir, "mr_all_result_12gene.RData"))
cat("\nMR主分析、异质性、多效性、留一法结果已保存\n")

#==================== 3.MR‑PRESSO：仅nSNP≥4的基因（CSE1L、EPB41L3），逐个基因运行 ====================
presso_list <- list()
presso_gene <- c("CSE1L","EPB41L3")
for(g in presso_gene){
  message("Running MR‑PRESSO for gene: ", g)
  dat_sub <- harmonised_fixed %>% filter(GeneSymbol == g)
  res_g <- mr_presso(
    BetaOutcome  = "beta.outcome",
    BetaExposure = "beta.exposure",
    SdOutcome    = "se.outcome",
    SdExposure   = "se.exposure",
    OUTLIERtest  = TRUE,
    DISTORTIONtest = TRUE,
    data         = dat_sub,
    NbDistribution = 1000
  )
  presso_list[[g]] <- res_g
}
save(presso_list, file = file.path(out_dir, "mr_presso_12gene.RData"))
print(presso_list)

#==================== 4.绘图：森林图、留一法图 ====================

cat("\n====所有MR分析完成，输出目录：",out_dir,"====\n")


#==================== 5.汇总：论文三线表（csv+Markdown） ====================
if(!requireNamespace("knitr", quietly = TRUE)) install.packages("knitr")
library(knitr)

ivw_tab <- mr_res_12 %>%
  filter(method == "Inverse variance weighted") %>%
  mutate(
    OR = exp(b),
    OR_l95 = exp(b - 1.96*se),
    OR_u95 = exp(b + 1.96*se),
    IVW_OR_CI = sprintf("%.3f (%.3f‑%.3f)", OR, OR_l95, OR_u95),
    IVW_P = signif(pval, 3)
  ) %>%
  select(GeneSymbol, n_snp_available, IVW_OR_CI, IVW_P)

# 异质性结果：提取IVW的Cochran Q，【删掉重复n_snp_available，防止join产生.x .y】
het_tab <- mr_het_12 %>%
  filter(method == "Inverse variance weighted") %>%
  mutate(Q_P = signif(Q_pval, 3)) %>%
  select(GeneSymbol, Q, Q_df, Q_P)   # 去掉 n_snp_available

# MR‑Egger截距多效性
egger_tab <- mr_pleio_12 %>%
  mutate(Egger_intercept_P = sprintf("%.3f (%.3g)", egger_intercept, pval)) %>%
  select(GeneSymbol, Egger_intercept_P)

# 主表合并：ivw_tab保留n_snp_available，het_tab不再带该列
tab_all <- ivw_tab %>%
  left_join(het_tab, by = "GeneSymbol") %>%
  left_join(egger_tab, by = "GeneSymbol")

# --------修复MR‑PRESSO取值--------
presso_summary <- tibble(
  GeneSymbol = names(presso_list),
  PRESSO_global_P = sapply(presso_list, function(x){
    p_val <- x$`MR-PRESSO results`$GlobalTest$Pvalue
    if(!is.null(p_val) && is.numeric(p_val) && !is.na(p_val)){
      return(signif(p_val, 3))
    }else{
      return(NA_real_)
    }
  })
)

tab_all <- tab_all %>%
  left_join(presso_summary, by = "GeneSymbol")

# 现在 n_snp_available 来自ivw_tab，列名正常，无.x/.y
tab_all <- tab_all %>%
  mutate(
    Egger_intercept_P = ifelse(n_snp_available < 3, "‑", Egger_intercept_P),
    PRESSO_global_P  = case_when(
      n_snp_available < 4 ~ "‑",
      is.na(PRESSO_global_P) ~ "‑",
      TRUE ~ as.character(PRESSO_global_P)
    ),
    Note = case_when(
      n_snp_available >=4 ~ "Full test available",
      n_snp_available ==3 ~ "Cannot perform MR‑PRESSO",
      n_snp_available ==2 ~ "Cannot perform MR‑Egger and MR‑PRESSO"
    )
  ) %>%
  arrange(desc(n_snp_available)) %>%
  rename(
    Gene = GeneSymbol,
    No_of_SNPs = n_snp_available,
    `IVW‑OR(95%CI)` = IVW_OR_CI,
    `IVW‑P` = IVW_P,
    `Cochran_Q_P` = Q_P,
    `MR‑Egger_intercept(P)` = Egger_intercept_P,
    `MR‑PRESSO_global_P` = PRESSO_global_P
  )

write.csv(tab_all, file.path(out_dir,"Table_MR_summary_12genes.csv"), row.names = FALSE)

cat("\n========论文三线表 Markdown格式，直接复制到文档========\n")
cat("```markdown\n")
print(kable(tab_all, format = "markdown"))
cat("```\n")
#####
library(dplyr)
library(ggplot2)

out_dir <- "D:/Bioinfor/colorectal/MR/eQTL"
mr_res_12 <- read.csv(file.path(out_dir, "mr_res_12gene.csv"))

# ---------------------- 1.提取IVW结果，计算OR与95%CI ----------------------
ivw_df <- mr_res_12 %>%
  filter(method == "Inverse variance weighted") %>%
  mutate(
    OR = exp(b),
    OR_l95 = exp(b - 1.96 * se),
    OR_u95 = exp(b + 1.96 * se),
    # 文本标签，图右侧展示
    or_text = sprintf("%.3f (%.3f‑%.3f)", OR, OR_l95, OR_u95),
    p_text = signif(pval,3)
  ) %>%
  arrange(OR) %>%   # 按OR升序排列
  mutate(GeneSymbol = factor(GeneSymbol, levels = unique(GeneSymbol)))

# ---------------------- 2.ggplot绘制森林图 ----------------------
p_forest_mr <- ggplot(ivw_df, aes(y = GeneSymbol, x = OR)) +
  # 无效线OR=1
  geom_vline(xintercept = 1, linetype = "dashed", colour = "#999999", linewidth = 0.8) +
  # 95%CI误差棒
  geom_errorbarh(aes(xmin = OR_l95, xmax = OR_u95), height = 0.3, colour = "#2c3e50", linewidth = 0.7) +
  # OR点估计
  geom_point(size = 3.2, colour = "#e74c3c") +
  # 右侧文字：OR(95%CI)
  geom_text(aes(label = or_text, x = max(ivw_df$OR_u95)*1.12),
            hjust = 0, size = 3.4) +
  labs(
    x = "Odds Ratio (OR)",
    y = "Gene symbol",
    title = "Mendelian randomization for 12 eQTL‑genes and colorectal cancer risk"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# 保存图片
ggsave(
  filename = file.path(out_dir, "MR_forest_IVW_ggplot.png"),
  plot = p_forest_mr,
  width = 14,
  height = 9,
  dpi = 300
)

print(p_forest_mr)
######
library(dplyr)
library(ggplot2)

out_dir <- "D:/Bioinfor/colorectal/MR/eQTL"
mr_loo <- read.csv(file.path(out_dir, "mr_leaveoneout_12gene.csv"))

# 只保留有留一法数据的基因（n>=3）
loo_df <- mr_loo %>%
  mutate(
    OR = exp(b),
    OR_l95 = exp(b - 1.96*se),
    OR_u95 = exp(b + 1.96*se),
    # 区分"All"完整估计和"剔除某SNP"
    label = ifelse(SNP == "All", "All SNPs", paste0("Drop ", SNP)),
    label = factor(label, levels = c("All SNPs", sort(unique(label[label != "All SNPs"]))))
  )

p_loo_facet <- ggplot(loo_df, aes(y = label, x = b)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey") +
  geom_errorbarh(aes(xmin = b - 1.96*se, xmax = b + 1.96*se),
                 height = 0.25, colour = "grey30") +
  geom_point(aes(colour = label == "All SNPs"), size = 2.5) +
  scale_colour_manual(values = c("FALSE" = "#2c3e50", "TRUE" = "#e74c3c"), guide = "none") +
  facet_wrap(~GeneSymbol, scales = "free_y", ncol = 2) +
  labs(x = "IVW beta (95% CI)", y = NULL,
       title = "Leave-one-out sensitivity analysis (genes with >=3 SNPs)") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_blank()
  )

ggsave(file.path(out_dir, "MR_leaveoneout_facet_ggplot.png"),
       plot = p_loo_facet, width = 14, height = 12, dpi = 300)
print(p_loo_facet)
#######------------漏斗图-----------------------##############################
library(dplyr)
library(ggplot2)
library(patchwork)

out_dir <- "D:/Bioinfor/colorectal/MR/eQTL"
load(file.path(out_dir, "harmonised_12gene.RData"))

harmonised_fixed <- harmonised_12gene %>%
  filter(
    id.outcome == "ebi-a-GCST90013862",
    !is.na(id.exposure),
    !is.na(id.outcome),
    mr_keep == TRUE
  ) %>%
  mutate(
    GeneSymbol = toupper(GeneSymbol),
    id.exposure = GeneSymbol,
    exposure = GeneSymbol,
    id.outcome = "ebi-a-GCST90013862",
    outcome = "Colorectal cancer"
  ) %>%
  distinct(GeneSymbol, SNP, .keep_all = TRUE)

funnel_gene <- harmonised_fixed %>%
  count(GeneSymbol) %>% filter(n >= 3) %>% pull(GeneSymbol)

# ============ 单基因漏斗图（重写，去掉ifelse嵌链和bquote） ============
plot_funnel_one <- function(dat_g, gene_name) {
  
  n_snp <- nrow(dat_g)
  
  # 计算Wald ratio + precision
  dat_funnel <- dat_g %>%
    mutate(
      wald_beta = beta.outcome / beta.exposure,
      wald_se   = se.outcome / abs(beta.exposure),
      precision = 1 / wald_se
    )
  
  # 调试：打印数据确认
  cat("\n=== ", gene_name, " (n=", n_snp, ") ===\n", sep="")
  print(dat_funnel %>% select(SNP, wald_beta, precision))
  
  # IVW估计
  ivw_res  <- TwoSampleMR::mr(dat_g, method_list = "mr_ivw")
  ivw_b    <- ivw_res$b[1]
  ivw_lo   <- ivw_b - 1.96 * ivw_res$se[1]
  ivw_hi   <- ivw_b + 1.96 * ivw_res$se[1]
  
  # Egger估计（单独tryCatch）
  egger_b <- tryCatch({
    e <- TwoSampleMR::mr(dat_g, method_list = "mr_egger_regression")
    e$b[1]
  }, error = function(e) NA)
  
  # 标签左右避让
  dat_funnel <- dat_funnel %>%
    mutate(
      label_x = ifelse(wald_beta < 0, wald_beta - 0.01, wald_beta + 0.01),
      hjust   = ifelse(wald_beta < 0, 1, 0)
    )
  
  # X轴范围：同时包含点和标签
  x_all <- c(dat_funnel$wald_beta, dat_funnel$label_x, ivw_lo, ivw_hi, 0)
  x_lim <- range(x_all, na.rm = TRUE)
  x_lim <- x_lim + c(-0.02, 0.02) * diff(x_lim)
  
  # 基础图层（不嵌ifelse）
  p <- ggplot(dat_funnel, aes(x = wald_beta, y = precision)) +
    annotate("rect",
             xmin = ivw_lo, xmax = ivw_hi,
             ymin = min(dat_funnel$precision) * 0.95,
             ymax = max(dat_funnel$precision) * 1.05,
             fill = "grey85", alpha = 0.6) +
    geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.4) +
    geom_vline(xintercept = ivw_b, colour = "#1f4e79", linewidth = 0.9) +
    geom_point(size = 3, shape = 21, fill = "white", colour = "black", stroke = 1) +
    geom_label(aes(x = label_x, label = SNP, hjust = hjust),
               vjust = 0.3, size = 2.8, colour = "black",
               label.size = 0, label.padding = unit(0.15, "lines")) +
    scale_x_continuous(limits = x_lim) +
    labs(
      title = paste0("Funnel: ", gene_name, " -> CRC (n=", n_snp, " SNPs)"),
      x = "Wald ratio (beta_outcome / beta_exposure)",
      y = "Precision (1 / SE_outcome)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title    = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.title    = element_text(size = 11),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3)
    )
  
  # Egger红线：单独加，不嵌在+链里
  if (!is.na(egger_b)) {
    p <- p + geom_vline(xintercept = egger_b, linetype = "dashed",
                        colour = "#c0392b", linewidth = 0.7)
  }
  
  return(p)
}

# ============ 循环绘图 ============
plot_list <- list()
for (g in funnel_gene) {
  cat("\n>>> Plotting ", g, "\n")
  dat_sub <- harmonised_fixed %>% filter(GeneSymbol == g)
  plot_list[[g]] <- plot_funnel_one(dat_sub, g)
}

# ============ 合并panel ============
panel_plot <- wrap_plots(plot_list, ncol = 2) +
  plot_annotation(
    title = "Funnel plots: eQTL genes -> colorectal cancer"
  ) &
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.margin = margin(t = 20, b = 10)
  )

# 保存
ggsave(file.path(out_dir, "MR_funnel_panel_5gene.png"),
       plot = panel_plot, width = 12, height = 15, dpi = 300)

print(panel_plot)
