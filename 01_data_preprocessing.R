rm(list = ls())
options(stringsAsFactors = FALSE, timeout = 3600)

work_dir <- "D:/Bioinfor/colorectal"
setwd(work_dir)
raw_dir <- file.path(work_dir, "GEO_CRC_raw")
result_dir <- file.path(work_dir, "results", "GSE196006_GSE146009_GSE110224")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

cran_packages <- c("BiocManager", "data.table", "ggplot2", "pheatmap")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

bioc_packages <- c(
  "GEOquery", "Biobase", "edgeR", "limma", "sva", "AnnotationDbi",
  "org.Hs.eg.db", "hgu133plus2.db"
)
missing_bioc <- bioc_packages[!vapply(
  bioc_packages, requireNamespace, logical(1), quietly = TRUE
)]
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(data.table)
  library(GEOquery)
  library(Biobase)
  library(edgeR)
  library(limma)
  library(sva)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(hgu133plus2.db)
  library(ggplot2)
  library(pheatmap)
})

# -----------------------------------------------------------------------------
# 1. GSE196006: RNA-seq raw counts, 21 paired normal/tumor samples
# -----------------------------------------------------------------------------
gse196_file <- file.path(raw_dir, "GSE196006_raw_counts.csv")


gse196 <- fread(gse196_file, check.names = FALSE)
gse196_id <- sub("\\..*$", "", as.character(gse196[[1]]))
gse196_counts <- as.matrix(gse196[, -1, with = FALSE])
storage.mode(gse196_counts) <- "numeric"
rownames(gse196_counts) <- gse196_id
colnames(gse196_counts) <- sub("_htseq\\.out$", "", colnames(gse196_counts))

group196 <- ifelse(
  grepl("[A-Z]0_", colnames(gse196_counts)), "Normal",
  ifelse(grepl("[A-Z]7_", colnames(gse196_counts)), "Tumor", NA_character_)
)
if (anyNA(group196)) stop("GSE196006 sample groups could not be completely identified")
#使用 edgeR 包创建 DGEList 对象
y196 <- edgeR::DGEList(counts = gse196_counts, group = factor(group196))

keep196 <- edgeR::filterByExpr(y196, group = factor(group196))
y196 <- y196[keep196, , keep.lib.sizes = FALSE]
y196 <- edgeR::calcNormFactors(y196, method = "TMM")
expr196 <- edgeR::cpm(y196, log = TRUE, prior.count = 1)

symbol196 <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = rownames(expr196),
  keytype = "ENSEMBL",
  column = "SYMBOL",
  multiVals = "first"
)
symbol196 <- toupper(trimws(as.character(symbol196)))
keep196_annotation <- !is.na(symbol196) & nzchar(symbol196)
expr196 <- expr196[keep196_annotation, , drop = FALSE]
symbol196 <- symbol196[keep196_annotation]
order196 <- unique(symbol196)
expr196_gene <- rowsum(expr196, group = symbol196, reorder = FALSE)
#统计`symbol196`中，按照`order196`规定顺序的各元素出现次数
count196 <- tabulate(match(symbol196, order196), nbins = length(order196))
#矩阵每一行的各个元素，除以该行对应的 `count196` 值
expr196_gene <- sweep(expr196_gene, 1, count196, "/")

info196 <- data.frame(
  sample_id = colnames(expr196_gene),
  dataset = "GSE196006",
  platform = "GPL15520",
  group = group196,
  stringsAsFactors = FALSE
)

rm(gse196, gse196_counts, y196, expr196)
invisible(gc())

# -----------------------------------------------------------------------------
# 2. GSE146009: RNA-seq raw counts from African American and Caucasian cohorts
# -----------------------------------------------------------------------------
gse146_aa_file <- file.path(raw_dir, "GSE146009_AfricanAmerican_RawCounts.txt")
gse146_ca_file <- file.path(raw_dir, "GSE146009_Caucasian_RawCounts.txt")


gse146_aa <- fread(gse146_aa_file, check.names = FALSE)
gse146_ca <- fread(gse146_ca_file, check.names = FALSE)
#将 gse146_aa 的第一列重命名为 GeneSymbol
setnames(gse146_aa, 1, "GeneSymbol")
setnames(gse146_ca, 1, "GeneSymbol")
gse146 <- merge(gse146_aa, gse146_ca, by = "GeneSymbol", all = FALSE)
symbol146 <- toupper(trimws(as.character(gse146$GeneSymbol)))
expr146_counts <- as.matrix(gse146[, -1, with = FALSE])
storage.mode(expr146_counts) <- "numeric"
rownames(expr146_counts) <- symbol146

group146 <- ifelse(
  grepl("T_S", colnames(expr146_counts), ignore.case = TRUE), "Tumor",
  ifelse(grepl("N_S", colnames(expr146_counts), ignore.case = TRUE),
         "Normal", NA_character_)
)
if (anyNA(group146)) stop("GSE146009 sample groups could not be completely identified")

valid146 <- !is.na(symbol146) & nzchar(symbol146) & !duplicated(symbol146)
expr146_counts <- expr146_counts[valid146, , drop = FALSE]
y146 <- edgeR::DGEList(counts = expr146_counts, group = factor(group146))
keep146 <- edgeR::filterByExpr(y146, group = factor(group146))
y146 <- y146[keep146, , keep.lib.sizes = FALSE]
y146 <- edgeR::calcNormFactors(y146, method = "TMM")
expr146_gene <- edgeR::cpm(y146, log = TRUE, prior.count = 1)

info146 <- data.frame(
  sample_id = colnames(expr146_gene),
  dataset = "GSE146009",
  platform = "GPL18573",
  group = group146,
  population = ifelse(grepl("^DS", colnames(expr146_gene)),
                      "African American", "Caucasian"),
  stringsAsFactors = FALSE
)

rm(gse146_aa, gse146_ca, gse146, expr146_counts, y146)
invisible(gc())

# -----------------------------------------------------------------------------
# 3. GSE110224: GPL570 Affymetrix microarray, 17 paired normal/tumor samples
# -----------------------------------------------------------------------------
gse110 <- GEOquery::getGEO(
  "GSE110224",
  GSEMatrix = TRUE,
  getGPL = FALSE,
  destdir = raw_dir
)
if (is.list(gse110)) gse110 <- gse110[[1]]
expr110 <- Biobase::exprs(gse110)
storage.mode(expr110) <- "numeric"
if (max(expr110, na.rm = TRUE) > 100) expr110 <- log2(expr110 + 1)
expr110 <- limma::normalizeBetweenArrays(expr110, method = "quantile")

probe110 <- rownames(expr110)
symbol110 <- AnnotationDbi::mapIds(
  hgu133plus2.db,
  keys = probe110,
  keytype = "PROBEID",
  column = "SYMBOL",
  multiVals = "first"
)
symbol110 <- toupper(trimws(as.character(symbol110)))
keep110 <- !is.na(symbol110) & nzchar(symbol110) &
  rowSums(is.finite(expr110)) == ncol(expr110)
expr110 <- expr110[keep110, , drop = FALSE]
symbol110 <- symbol110[keep110]
order110 <- unique(symbol110)
expr110_gene <- rowsum(expr110, group = symbol110, reorder = FALSE)
count110 <- tabulate(match(symbol110, order110), nbins = length(order110))
expr110_gene <- sweep(expr110_gene, 1, count110, "/")
#从 GSE110 的 ExpressionSet 对象中，获取样本的表型 / 临床信息
pd110 <- Biobase::pData(gse110)
title110 <- as.character(pd110[colnames(expr110_gene), "title"])
group110 <- ifelse(
  grepl("cancer", title110, ignore.case = TRUE), "Tumor",
  ifelse(grepl("norm", title110, ignore.case = TRUE), "Normal", NA_character_)
)
if (anyNA(group110)) stop("GSE110224 sample groups could not be completely identified")

info110 <- data.frame(
  sample_id = colnames(expr110_gene),
  sample_title = title110,
  dataset = "GSE110224",
  platform = "GPL570",
  group = group110,
  stringsAsFactors = FALSE
)

rm(gse110, expr110)
invisible(gc())

# -----------------------------------------------------------------------------
# 4. Common genes, cross-platform normalization, merging and ComBat correction
# -----------------------------------------------------------------------------
common_genes <- Reduce(
  intersect,
  list(rownames(expr196_gene), rownames(expr146_gene), rownames(expr110_gene))
)
if (length(common_genes) < 1000) {
  stop("Only ", length(common_genes), " common genes were found")
}

expr196_gene <- expr196_gene[common_genes, , drop = FALSE]
expr146_gene <- expr146_gene[common_genes, , drop = FALSE]
expr110_gene <- expr110_gene[common_genes, , drop = FALSE]
merged <- cbind(expr196_gene, expr146_gene, expr110_gene)

sample_info <- rbind(
  info196[, c("sample_id", "dataset", "platform", "group")],
  info146[, c("sample_id", "dataset", "platform", "group")],
  info110[, c("sample_id", "dataset", "platform", "group")]
)
sample_info <- sample_info[match(colnames(merged), sample_info$sample_id), ]
if (anyNA(sample_info$sample_id)) stop("Sample information does not match merged matrix")

# Quantile normalization puts RNA-seq log-CPM and microarray log2 expression on
# a comparable empirical distribution before batch correction.
merged_normalized <- limma::normalizeBetweenArrays(merged, method = "quantile")
#对归一化表达矩阵，逐行计算每个基因的表达方差
gene_variance <- apply(merged_normalized, 1, var)
merged_normalized <- merged_normalized[
  is.finite(gene_variance) & gene_variance > 0, , drop = FALSE
]
#把 dataset 转为因子，手动设定 3 个 GSE 数据集的水平顺序。
sample_info$dataset <- factor(
  sample_info$dataset,
  levels = c("GSE196006", "GSE146009", "GSE110224")
)
#将sample_info数据框的group列转为因子类型
sample_info$group <- factor(sample_info$group, levels = c("Normal", "Tumor"))

combat_design <- model.matrix(~ group, data = sample_info)
merged_combat <- sva::ComBat(
  dat = merged_normalized,
  batch = sample_info$dataset,
  mod = combat_design,
  par.prior = TRUE,
  prior.plots = FALSE
)

write.csv(
  data.frame(GeneSymbol = rownames(merged_normalized), merged_normalized,
             check.names = FALSE),
  file.path(result_dir, "merged_before_ComBat.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(GeneSymbol = rownames(merged_combat), merged_combat,
             check.names = FALSE),
  file.path(result_dir, "merged_after_ComBat.csv"),
  row.names = FALSE
)
write.csv(sample_info, file.path(result_dir, "sample_information.csv"),
          row.names = FALSE)

# -----------------------------------------------------------------------------
# 5. PCA before and after batch correction
# 5.1 读回 write.csv 生成的三个文件
# 第1列是 GeneSymbol，作为行名；check.names=FALSE 保留样本名的特殊字符
expr_before <- read.csv(file.path(result_dir, "merged_before_ComBat.csv"),
                        row.names = 1, check.names = FALSE)
expr_after  <- read.csv(file.path(result_dir, "merged_after_ComBat.csv"),
                        row.names = 1, check.names = FALSE)
sample_info <- read.csv(file.path(result_dir, "sample_information.csv"),
                        stringsAsFactors = FALSE)

# 检查样本名与 sample_info 顺序一致（gene 行数、sample 列数）
ord <- match(sample_info$sample, colnames(expr_before))
expr_before <- expr_before[, ord]
expr_after  <- expr_after[, ord]
stopifnot(ncol(expr_before) == nrow(sample_info),
          ncol(expr_after) == nrow(sample_info))
# -----------------------------------------------------------------------------
# 5.2 高变异基因（前 5000，与之前一致）
top_n <- min(5000, nrow(expr_before))
top_genes <- names(sort(apply(expr_before, 1, var), decreasing = TRUE))[
  seq_len(top_n)
]

# 5.3 两个阶段分别做 PCA
pca_before <- prcomp(t(expr_before[top_genes, , drop = FALSE]),
                     center = TRUE, scale. = FALSE)
pca_after  <- prcomp(t(expr_after[top_genes, , drop = FALSE]),
                     center = TRUE, scale. = FALSE)

var_before <- 100 * pca_before$sdev^2 / sum(pca_before$sdev^2)
var_after  <- 100 * pca_after$sdev^2 / sum(pca_after$sdev^2)

df_pca_before <- cbind(sample_info, PC1 = pca_before$x[,1], PC2 = pca_before$x[,2])
df_pca_after  <- cbind(sample_info, PC1 = pca_after$x[,1],  PC2 = pca_after$x[,2])

pca_data_before <- cbind(
  sample_info,
  PC1 = pca_before$x[, 1],
  PC2 = pca_before$x[, 2],
  stage = factor("Before batch correction",
                 levels = c("Before batch correction", "After batch correction")),
  PC1_var = var_before[1],
  PC2_var = var_before[2]
)
pca_data_after <- cbind(
  sample_info,
  PC1 = pca_after$x[, 1],
  PC2 = pca_after$x[, 2],
  stage = factor("After batch correction",
                 levels = c("Before batch correction", "After batch correction")),
  PC1_var = var_after[1],
  PC2_var = var_after[2]
)
pca_data <- rbind(pca_data_before, pca_data_after)

write.csv(pca_data, file.path(result_dir, "PCA_coordinates.csv"),
          row.names = FALSE)
# ===================== 绘图（文献风格） =====================

#----------- A: PCA Before batch correction -----------
pca_before_plot <- ggplot(df_pca_before, aes(x = PC1, y = PC2)) +
  stat_ellipse(aes(color = dataset, fill = dataset), geom = "polygon", alpha = 0.2, level = 0.95) +
  geom_point(aes(color = dataset, shape = dataset), size = 2.2, alpha = 0.9) +
  labs(title = "Before batch correction", x = "PC1", y = "PC2") +
  scale_color_manual(name = "Type",
                     values = c(GSE149507="#E74C3C", GSE43346="#27AE60", GSE73160="#3498DB")) +
  scale_fill_manual(name = "Type",
                    values = c(GSE149507="#E74C3C", GSE43346="#27AE60", GSE73160="#3498DB")) +
  scale_shape_manual(name = "Type",
                     values = c(GSE149507=16, GSE43346=17, GSE73160=15)) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color="grey90", linewidth=0.2),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size =12))

#----------- B: PCA After batch correction -----------
pca_after_plot <- ggplot(df_pca_after, aes(x = PC1, y = PC2)) +
  stat_ellipse(aes(color = dataset, fill = dataset), geom = "polygon", alpha = 0.2, level = 0.95) +
  geom_point(aes(color = dataset, shape = dataset), size = 2.2, alpha = 0.9) +
  labs(title = "After batch correction", x = "PC1", y = "PC2") +
  scale_color_manual(name = "Type",
                     values = c(GSE149507="#E74C3C", GSE43346="#27AE60", GSE73160="#3498DB")) +
  scale_fill_manual(name = "Type",
                    values = c(GSE149507="#E74C3C", GSE43346="#27AE60", GSE73160="#3498DB")) +
  scale_shape_manual(name = "Type",
                     values = c(GSE149507=16, GSE43346=17, GSE73160=15)) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color="grey90", linewidth=0.2),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size =12))

# -----------------------------------------------------------------------------
# 6. Differential expression, volcano plot and heatmap
# -----------------------------------------------------------------------------
# Tumor versus Normal is tested on the ComBat-corrected common-gene matrix.
# The design retains the biological group while the dataset batch has already
# been corrected above.
de_design <- model.matrix(~ 0 + sample_info$group)
colnames(de_design) <- c("Normal", "Tumor")
de_contrast <- limma::makeContrasts(Tumor - Normal, levels = de_design)
de_fit <- limma::lmFit(merged_combat, de_design)
de_fit <- limma::contrasts.fit(de_fit, de_contrast)
de_fit <- limma::eBayes(de_fit, robust = TRUE)
de_table <- limma::topTable(
  de_fit,
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)
de_table$GeneSymbol <- rownames(de_table)
de_table$negLog10FDR <- -log10(pmax(de_table$adj.P.Val, .Machine$double.xmin))
de_table$Category <- "Not significant"
de_table$Category[de_table$adj.P.Val < 0.05 & de_table$logFC >= 1] <- "Up"
de_table$Category[de_table$adj.P.Val < 0.05 & de_table$logFC <= -1] <- "Down"
de_table <- de_table[, c(
  "GeneSymbol", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val",
  "B", "negLog10FDR", "Category"
)]
write.csv(
  de_table,
  file.path(result_dir, "DEG_Tumor_vs_Normal_limma.csv"),
  row.names = FALSE
)
write.csv(
  de_table[de_table$Category != "Not significant", , drop = FALSE],
  file.path(result_dir, "DEG_Tumor_vs_Normal_significant.csv"),
  row.names = FALSE
)

# ===== 从已保存的 CSV 直接读取 de_table =====
de_table <- read.csv(file.path(result_dir, "DEG_Tumor_vs_Normal_limma.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)

# 类型安全：把 Category 设为因子，保证图例顺序稳定、drop=FALSE 生效
de_table$Category <- factor(de_table$Category,
                            levels = c("Up", "Down", "Not significant"))

de_table$Category <- factor(de_table$Category,
                            levels = c("Up", "Down", "Not significant"))

volcano_plot <- ggplot(
  de_table,
  aes(x = logFC, y = negLog10FDR, color = Category)
) +
  geom_point(size = 1.4, alpha = 0.72, na.rm = TRUE) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey45") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey45") +
  scale_color_manual(
    values = c(
      Down = "#0072B2",
      `Not significant` = "grey72",
      Up = "#D55E00"
    ),
    drop = FALSE
  ) +
  labs(
    title = "GSE196006 + GSE146009 + GSE110224",
    subtitle = "Tumor versus Normal after ComBat correction",
    x = "log2 fold change",
    y = "-log10 adjusted P value",
    color = "Category"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom")

ggsave(
  file.path(result_dir, "Volcano_Tumor_vs_Normal.png"), volcano_plot,
  width = 7.5, height = 6.2, dpi = 300
)
ggsave(
  file.path(result_dir, "Volcano_Tumor_vs_Normal.pdf"), volcano_plot,
  width = 7.5, height = 6.2
)

# 输出统计信息便于核对
cat("Total DEG rows:", nrow(de_table), "\n")
print(table(de_table$Category))


sig_de <- de_table[de_table$Category != "Not significant", , drop = FALSE]
sig_de <- sig_de[order(sig_de$adj.P.Val, -abs(sig_de$logFC)), , drop = FALSE]
if (nrow(sig_de) > 0) {
  up_genes <- sig_de$GeneSymbol[sig_de$Category == "Up"]
  down_genes <- sig_de$GeneSymbol[sig_de$Category == "Down"]
  up_genes <- head(up_genes, 50)
  down_genes <- head(down_genes, 50)
  heat_genes <- unique(c(up_genes, down_genes))
} else {
  heat_genes <- rownames(de_table)[order(de_table$adj.P.Val)]
  heat_genes <- head(heat_genes, min(50, length(heat_genes)))
}
#heat_genes <- heat_genes[heat_genes %in% rownames(expr_after)]
heat_genes <- head(heat_genes, 40)
heat_matrix <- expr_after[heat_genes, , drop = FALSE]
heat_matrix <- t(scale(t(heat_matrix)))
heat_matrix[!is.finite(heat_matrix)] <- 0
colnames(heat_matrix) <- sample_info$sample_id
annotation_col <- data.frame(
  Group = sample_info$group,
  Dataset = sample_info$dataset,
  row.names = sample_info$sample_id,
  check.names = FALSE
)
annotation_colors <- list(
  Group = c(Normal = "#0072B2", Tumor = "#D55E00"),
  Dataset = c(
    GSE196006 = "#009E73",
    GSE146009 = "#CC79A7",
    GSE110224 = "#E69F00"
  )
)
write.csv(
  data.frame(GeneSymbol = rownames(heat_matrix), heat_matrix,
             check.names = FALSE),
  file.path(result_dir, "DEG_heatmap_Zscore_matrix.csv"),
  row.names = FALSE
)

# ========== pheatmap绘图，同时转为ggplot对象用于patchwork拼图 ==========
# silent=TRUE：不在设备直接绘图，捕获绘图对象
pheat_obj <- pheatmap::pheatmap(
  heat_matrix,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 5,
  border_color = NA,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  silent = TRUE
)

# 转ggplot对象，用于patchwork拼图（D子图）
library(pheatmap)
library(ggplotify)
heatmap_plot <- as.ggplot(pheat_obj)

# 输出独立png文件
pheatmap::pheatmap(
  heat_matrix,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 6,
  border_color = NA,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  filename = file.path(result_dir, "Heatmap_DEG_Tumor_vs_Normal.png"),
  width = 12,
  height = max(7, 0.16 * nrow(heat_matrix) + 3),
  res = 300
)
# 输出独立pdf文件
pheatmap::pheatmap(
  heat_matrix,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 6,
  border_color = NA,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  filename = file.path(result_dir, "Heatmap_DEG_Tumor_vs_Normal.pdf"),
  width = 12,
  height = max(7, 0.16 * nrow(heat_matrix) + 3)
)

saveRDS(
  list(
    DEG = de_table,
    significant_DEG = sig_de,
    heatmap_matrix = heat_matrix,
    heatmap_annotation = annotation_col
  ),
  file.path(result_dir, "DEG_Tumor_vs_Normal_results.rds"),
  compress = "xz"
)
cat("DEGs with FDR < 0.05 and |logFC| >= 1:", nrow(sig_de), "\n")

####拼图
library(patchwork)
figure_all <- (pca_before_plot + pca_after_plot) /
  (volcano_plot + heatmap_plot) +
  plot_annotation(
    tag_levels = "A",
    tag_prefix = "(",
    tag_suffix = ")"
  ) +
  plot_layout(
    nrow = 2,
    guides = "collect",
    widths = c(1, 1)
  ) &
  theme(plot.tag = element_text(size = 13, face = "bold"))

ggsave(
  "Full_Figure_ABCD_heatmap_revised.png",
  figure_all,
  width = 14,
  height = 10,
  dpi = 300
)

