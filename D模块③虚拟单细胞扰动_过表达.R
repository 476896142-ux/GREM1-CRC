# ================================================================
# scTenifold 过表达分析（GREM1 OE）—— 高变基因加速版
# 原理：WT 网络 vs GREM1×10 过表达网络，流形对齐找响应基因
# ================================================================

library(scTenifoldNet)
library(scTenifoldKnk)
library(dplyr)
library(Seurat)

# ===== 目录设置 =====
setwd("D:/Bioinfor/colorectal")
base_dir <- "D:/Bioinfor/colorectal/single_raw"
raw_dir <- file.path(base_dir, "raw")
result_dir <- file.path(base_dir, "results")
figure_dir <- file.path(result_dir, "figures")

# ================================================================
# 1. 读取 Seurat 对象 + 细胞过滤
# ================================================================
sc <- readRDS(file.path(result_dir, "seurat_merged.rds"))

# 去掉 celltype 为 NA 或空的细胞
sc <- subset(sc, cells = colnames(sc)[!is.na(sc$celltype) & sc$celltype != ""])
sc$celltype <- droplevels(as.factor(sc$celltype))
cat("过滤后细胞数:", ncol(sc), "\n")
cat("细胞类型:", paste(levels(sc$celltype), collapse = ", "), "\n")

# ================================================================
# 2. 提取基质细胞计数矩阵 + 高变基因筛选
# ================================================================
# 从过滤后的 sc 提取原始计数（用 count_matrix 避免和 R 内置 counts() 冲突）
count_matrix <- GetAssayData(sc, layer = "counts")
cat("计数矩阵:", nrow(count_matrix), "基因 ×", ncol(count_matrix), "细胞\n")

# 提取基质细胞（GREM1 主要表达的细胞类型）
stromal_cells <- colnames(sc)[sc$celltype == "Stromal cells"]
counts_wt <- count_matrix[, stromal_cells]
cat("基质细胞:", ncol(counts_wt), "个\n")

# 只保留高变基因 + GREM1（和敲除分析一致，加速计算）
variable_genes <- VariableFeatures(sc)
variable_genes <- unique(c(variable_genes, "GREM1"))
counts_wt_fast <- counts_wt[rownames(counts_wt) %in% variable_genes, ]
cat("加速版基因数:", nrow(counts_wt_fast), "\n")

# ================================================================
# 3. 构建 WT（野生型）基因调控网络
# ================================================================
cat("\n=== 构建 WT 网络 ===\n")
cat("开始时间:", as.character(Sys.time()), "\n")
set.seed(123)

args(scTenifoldNet::makeNetworks)

# ===== 构建 WT 网络 =====
wt_net <- scTenifoldNet::makeNetworks(
  X = counts_wt_fast,        # 修正：X 不是 countMatrix
  nNet = 3,
  nCells = 500,
  scaleScores = TRUE,
  symmetric = FALSE,
  q = 0.95,                  # 默认 0.95
  nCores = 16                # 你的 CPU 核数
)

#对单细胞张量网络 wt_net 执行张量分解
wt_tensor <- scTenifoldNet::tensorDecomposition(
  wt_net, maxIter = 1000, maxError = 1e-5
)
cat("WT 网络完成:", as.character(Sys.time()), "\n")

# ================================================================

# ===== 模拟过表达 =====
counts_oe <- counts_wt_fast
#定位表达矩阵中 GREM1 基因的行位置，用于后续提取该基因的表达数据
grem1_idx <- which(rownames(counts_oe) == "GREM1")
counts_oe[grem1_idx, ] <- counts_oe[grem1_idx, ] * 10

# ===== 构建 OE 网络 =====
oe_net <- scTenifoldNet::makeNetworks(
  X = counts_oe,
  nNet = 3,
  nCells = 500,
  scaleScores = TRUE,
  symmetric = FALSE,
  q = 0.95,
  nCores = 16
)
oe_tensor <- scTenifoldNet::tensorDecomposition(oe_net, maxIter = 1000, maxError = 1e-5)

# ===== 流形对齐 =====

ma_result <- scTenifoldNet::manifoldAlignment(
  X = wt_tensor$X,
  Y = oe_tensor$X,
  d = 2,
  nCores = 16
)

# 7. 提取响应基因（清理基因名 + 去重）
# ================================================================
# 计算欧氏距离，去掉 X_/Y_ 前缀，重复基因取距离大的
gene_dist <- data.frame(
  gene_raw = rownames(ma_result),
  gene = sub("^[XY]_", "", rownames(ma_result)),   # 去掉 X_ 和 Y_ 前缀
  NLMA1 = ma_result[, 1],
  NLMA2 = ma_result[, 2],
  distance = sqrt(ma_result[, 1]^2 + ma_result[, 2]^2)
) %>%
  group_by(gene) %>%
  summarise(distance = max(distance), .groups = "drop") %>%
  arrange(desc(distance))

cat("清理后基因数:", nrow(gene_dist), "\n")

# 排除 GREM1 本身，取 TOP20
top20_oe <- gene_dist %>%
  filter(gene != "GREM1") %>%
  head(20)

cat("\n=== GREM1 过表达 TOP20 响应基因 ===\n")
print(as.data.frame(top20_oe), row.names = FALSE)

# ================================================================
# 8. 和敲除结果对比
# ================================================================
knk_top20 <- read.csv(file.path(result_dir, "GREM1_TOP20_response_genes.csv"))

cat("\n=== 敲除 vs 过表达 对比 ===\n")
cat("敲除 TOP10:", paste(head(knk_top20$gene, 10), collapse = ", "), "\n")
cat("过表达 TOP10:", paste(head(top20_oe$gene, 10), collapse = ", "), "\n")

common_genes <- intersect(knk_top20$gene, top20_oe$gene)
cat("共同基因数:", length(common_genes), "\n")
if (length(common_genes) > 0) {
  cat("共同基因:", paste(common_genes, collapse = ", "), "\n")
}

# ================================================================
# 9. 保存结果
# ================================================================
write.csv(gene_dist,
          file.path(result_dir, "GREM1_Overexpression_all_genes.csv"),
          row.names = FALSE)
write.csv(top20_oe,
          file.path(result_dir, "GREM1_Overexpression_TOP20.csv"),
          row.names = FALSE)
saveRDS(list(wt = wt_tensor, oe = oe_tensor, ma = ma_result, genes = gene_dist),
        file.path(result_dir, "scTenifoldOe_GREM1.rds"))

cat("\n========================================\n")
cat("过表达分析完成！结果保存在:", result_dir, "\n")
cat("  - GREM1_Overexpression_all_genes.csv\n")
cat("  - GREM1_Overexpression_TOP20.csv\n")
cat("  - scTenifoldOe_GREM1.rds\n")
cat("========================================\n")
########################################################################
library(clusterProfiler)
library(org.Hs.eg.db)

# ===== 用 TOP20 过表达基因做 ORA 富集 =====
top20_oe <- read.csv(file.path(result_dir, "GREM1_Overexpression_TOP20.csv"))

# ID 转换
eg_oe <- bitr(top20_oe$gene, fromType = "SYMBOL",
              toType = "ENTREZID", OrgDb = org.Hs.eg.db)
cat("映射成功:", nrow(eg_oe), "/ 20\n")

# GO 富集
go_oe <- enrichGO(gene = eg_oe$ENTREZID, OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID", ont = "ALL",
                  pvalueCutoff = 1, qvalueCutoff = 1)   # 阈值放宽到1

go_oe_df <- as.data.frame(go_oe) %>%
  arrange(pvalue) %>%
  head(15)

cat("\n=== 过表达 TOP20 GO 富集 TOP15 ===\n")
print(go_oe_df[, c("Description","ONTOLOGY","pvalue","p.adjust")], row.names = FALSE)

# KEGG 富集
kegg_oe <- enrichKEGG(gene = eg_oe$ENTREZID,
                      pvalueCutoff = 1, qvalueCutoff = 1)

kegg_oe_df <- as.data.frame(kegg_oe) %>%
  arrange(pvalue) %>%
  head(15)

cat("\n=== 过表达 TOP20 KEGG 富集 TOP15 ===\n")
print(kegg_oe_df[, c("Description","pvalue","p.adjust")], row.names = FALSE)

# 保存
write.csv(as.data.frame(go_oe), file.path(result_dir, "GREM1_OE_GO_ORA.csv"), row.names = FALSE)
write.csv(as.data.frame(kegg_oe), file.path(result_dir, "GREM1_OE_KEGG_ORA.csv"), row.names = FALSE)

# 画图（如果有通路）
if (nrow(go_oe_df) > 0) {
  p <- dotplot(go_oe, showCategory = 10) +
    facet_grid(ONTOLOGY ~ ., scales = "free_y") +
    ggtitle("GREM1 Overexpression TOP20 - GO Enrichment")
  ggsave(file.path(figure_dir, "12_GREM1_OE_GO.png"), p, width = 10, height = 8)
}
if (nrow(kegg_oe_df) > 0) {
  p2 <- dotplot(kegg_oe, showCategory = 10) +
    ggtitle("GREM1 Overexpression TOP20 - KEGG")
  ggsave(file.path(figure_dir, "13_GREM1_OE_KEGG.png"), p2, width = 8, height = 6)
}
cat("\n富集分析完成！\n")
library(ggplot2)
library(ggrepel)  # 用于基因名标注，避免重叠
##########################火山图
