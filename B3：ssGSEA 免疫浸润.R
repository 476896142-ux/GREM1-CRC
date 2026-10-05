# ===== 安装包（首次运行） =====
BiocManager::install("GSVA")
BiocManager::install("GSEABase")
library(GSVA)
library(GSEABase)
library(dplyr)
library(tibble)
library(pheatmap)

# ===== 1. 免疫细胞 marker 基因集（16 种主要免疫细胞） =====
immune_sets <- list(
  `B cells`         = c("MS4A1","CD19","CD79A","CD79B","BLK","BANK1"),
  `CD8+ T cells`    = c("CD8A","CD8B","GZMA","GZMB","PRF1","CD247"),
  `CD4+ T cells`    = c("CD4","IL7R","LTB","TRAT1","CD3D","CD3E"),
  `Tregs`           = c("FOXP3","IL2RA","CTLA4","TIGIT","IKZF2"),
  `Th1 cells`       = c("TBX21","IFNG","STAT1","CXCR3","TNF"),
  `Th2 cells`       = c("GATA3","IL4","IL5","IL13","STAT6"),
  `Th17 cells`      = c("RORC","IL17A","IL17F","CCR6","IL21"),
  `NK cells`        = c("KLRD1","NKG7","GNLY","KLRF1","NCR1"),
  `Macrophages`     = c("CD68","CD163","MRC1","MSR1","CSF1R"),
  `M1 Macrophages`  = c("NOS2","IL1B","TNF","CXCL9","CXCL10"),
  `M2 Macrophages`  = c("CD163","MRC1","MSR1","ARG1","IL10"),
  `Neutrophils`     = c("FCGR3B","S100A8","S100A9","CEACAM8","CXCR2"),
  `Dendritic cells` = c("CD1C","THBD","CLEC9A","XCR1","FCER1A"),
  `Mast cells`      = c("TPSAB1","TPSB2","CPA3","KIT","HDC"),
  `Monocytes`       = c("CD14","FCGR3A","S100A8","S100A9","LYZ"),
  `Eosinophils`     = c("CLC","IL5RA","PRG2","SIGLEC8","RNASE3")
)

# 转成 GeneSetCollection
immune_gsc <- GeneSetCollection(mapply(function(genes, name) {
  GeneSet(geneIds = genes, setName = name)
}, immune_sets, names(immune_sets)))

# ===== 2. ssGSEA 计算（仅肿瘤样本） =====
tumor_samples <- subtype$sample   # 376 个肿瘤样本
expr_tumor <- expr[, tumor_samples, drop = FALSE]   # 基因 × 样本

# 新版 API：先用 ssgseaParam 构建参数，geneSets 直接用命名 list
param <- ssgseaParam(
  expr     = as.matrix(expr_tumor),
  geneSets = immune_sets    # 直接用上面定义的 list，不需要 GeneSetCollection
)
set.seed(123)
ssgsea_scores <- gsva(param, verbose = FALSE)

# 结果：16 种免疫细胞 × 376 样本
dim(ssgsea_scores)
head(ssgsea_scores[, 1:3])

# ===== 3. 合并亚型标签 =====
ssgsea_df <- as.data.frame(t(ssgsea_scores)) %>%
  rownames_to_column("sample") %>%
  left_join(subtype, by = "sample")

# ===== 4. 亚型间免疫浸润差异（Kruskal-Wallis 检验） =====
immune_cells <- names(immune_sets)
kw_results <- data.frame(
  ImmuneCell = immune_cells,
  p_value = sapply(immune_cells, function(cell) {
    kruskal.test(as.formula(paste0("`", cell, "` ~ Subtype")), data = ssgsea_df)$p.value
  })
) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

print(kw_results)
cat("\n显著差异的免疫细胞数（FDR<0.05）:", sum(kw_results$p_adj < 0.05), "\n")

# ===== 5. 热图：各亚型免疫浸润谱 =====
# 按亚型排序
ssgsea_ordered <- ssgsea_df %>% arrange(Subtype)
mat <- as.matrix(ssgsea_ordered[, immune_cells])
rownames(mat) <- ssgsea_ordered$sample

# 亚型注释条
anno_col <- data.frame(Subtype = factor(ssgsea_ordered$Subtype))
rownames(anno_col) <- ssgsea_ordered$sample

png("MR/results/ssgsea_immune_heatmap.png",
    width = 2000, height = 1200, res = 200)
pheatmap(t(mat),
         annotation_col = anno_col,
         show_colnames = FALSE,
         scale = "row",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         main = "ssGSEA immune infiltration across 4 subtypes")
dev.off()

# ===== 6. 保存结果 =====
write.csv(ssgsea_df, "MR/results/ssgsea_scores_by_subtype.csv", row.names = FALSE)


write.csv(kw_results, "MR/results/ssgsea_kw_test.csv", row.names = FALSE)



cat("已保存: ssgsea_scores_by_subtype.csv, ssgsea_kw_test.csv, ssgsea_immune_heatmap.png\n")


# 各亚型的平均免疫浸润分数
subtype_immune_mean <- ssgsea_df %>%
  group_by(Subtype) %>%
  summarise(across(all_of(immune_cells), mean, na.rm = TRUE),
            n = n(), .groups = "drop")

print(as.data.frame(subtype_immune_mean), row.names = FALSE)

# 画 top 6 免疫细胞的箱线图
library(ggplot2)
library(tidyr)

top6 <- kw_results$ImmuneCell[1:6]
ssgsea_long <- ssgsea_df %>%
  dplyr::select(sample, Subtype, dplyr::all_of(top6)) %>%
  tidyr::pivot_longer(cols = dplyr::all_of(top6),
                      names_to = "ImmuneCell", values_to = "Score")

p <- ggplot(ssgsea_long, aes(x = factor(Subtype), y = Score, fill = factor(Subtype))) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7) +
  facet_wrap(~ ImmuneCell, scales = "free_y", ncol = 3) +
  labs(x = "Subtype", y = "ssGSEA score", title = "Top 6 immune cells across 4 subtypes") +
  theme_bw() +
  theme(legend.position = "none",
        strip.text = element_text(size = 10),
        axis.text.x = element_text(size = 9))

ggsave("MR/results/ssgsea_boxplot_top6.png", p, width = 12, height = 8, dpi = 200)
cat("已保存: ssgsea_boxplot_top6.png\n")

