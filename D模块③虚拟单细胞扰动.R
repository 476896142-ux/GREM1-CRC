# =============================================
# 路径设置：按用户指定目录组织原始文件、结果和图形。
setwd("D:/Bioinfor/colorectal")
base_dir <- "D:/Bioinfor/colorectal/single_raw"
raw_dir <- file.path(base_dir, "raw")
result_dir <- file.path(base_dir, "results")
figure_dir <- file.path(result_dir, "figures")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

library(scTenifoldKnk)
library(Seurat)
library(dplyr)

sc <- readRDS(file.path(result_dir, "seurat_merged.rds"))
sc <- subset(sc, cells = colnames(sc)[!is.na(sc$celltype) & sc$celltype != ""])
sc$celltype <- droplevels(as.factor(sc$celltype))

# 提取基质细胞计数矩阵
stromal_cells <- colnames(sc)[sc$celltype == "Stromal cells"]
counts <- GetAssayData(sc, layer = "counts")
counts_stromal <- counts[, stromal_cells]
cat("基质细胞:", ncol(counts_stromal), "个\n")

# 只保留高变基因 + GREM1（加速）
variable_genes <- VariableFeatures(sc)
variable_genes <- unique(c(variable_genes, "GREM1"))
counts_fast <- counts_stromal[rownames(counts_stromal) %in% variable_genes, ]
cat("加速版基因数:", nrow(counts_fast), "\n")

# 看 CPU 核数
n_cores <- parallel::detectCores()
cat("CPU 核数:", n_cores, "\n")

# 跑 scTenifoldKnk
set.seed(123)
cat("开始时间:", as.character(Sys.time()), "\n")

knk_result <- scTenifoldKnk(
  countMatrix = counts_fast,
  gKO = "GREM1",
  qc_minPCT = 0.01,
  nc_nNet = 3,
  nc_lambda = 0.1,
  nCores = n_cores    # 自动用全部核
)

cat("结束时间:", as.character(Sys.time()), "\n")
# ===== 1. 提取全部差异结果，按 distance 排序 =====
all_genes <- knk_result$diffRegulation %>%
  arrange(desc(distance))

cat("全部基因数:", nrow(all_genes), "\n")
cat("显著基因数 (p.adj<0.05):", sum(all_genes$p.adj < 0.05), "\n")

# ===== 2. 取 TOP20 响应基因（排除 GREM1 本身）=====
top20 <- all_genes %>%
  filter(gene != "GREM1") %>%
  head(20)

cat("\n=== GREM1 敲除 TOP20 响应基因 ===\n")
print(as.data.frame(top20[, c("gene","distance","Z","FC","p.value")]), row.names = FALSE)

write.csv(top20, file.path(result_dir, "GREM1_TOP20_response_genes.csv"), row.names = FALSE)

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
# 全部基因按 distance 排序（预排序 GSEA）
all_genes_sorted <- knk_result$diffRegulation %>%
  arrange(desc(distance)) %>%
  filter(gene != "GREM1")

# GeneSymbol → EntrezID
gene_map <- bitr(all_genes_sorted$gene, fromType = "SYMBOL",
                 toType = "ENTREZID", OrgDb = org.Hs.eg.db)
all_genes_sorted <- all_genes_sorted %>%
  inner_join(gene_map, by = c("gene" = "SYMBOL"))

# 构建预排序基因列表
gene_list <- all_genes_sorted$distance
names(gene_list) <- all_genes_sorted$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

# GO GSEA
gse_go <- gseGO(geneList = gene_list, OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID", ont = "BP",
                pvalueCutoff = 0.5, verbose = FALSE)

# 筛选 p.adjust < 0.25，按 NES 排序
gse_go_sig <- as.data.frame(gse_go) %>%
  filter(p.adjust < 0.25) %>%
  arrange(desc(NES))

cat("GO 显著通路数 (p.adj<0.25):", nrow(gse_go_sig), "\n")

if (nrow(gse_go_sig) > 0) {
  cat("\n--- TOP15 正向富集通路 ---\n")
  print(head(gse_go_sig[, c("Description","NES","pvalue","p.adjust")], 15),
        row.names = FALSE)
  
  # ECM/胶原相关通路
  ecm_go <- gse_go_sig[grep("collagen|extracellular matrix|ECM|matrix|fibroblast|connective tissue|skin|epiderm",
                            gse_go_sig$Description, ignore.case = TRUE), ]
  if (nrow(ecm_go) > 0) {
    cat("\n--- ECM/胶原相关通路 ---\n")
    print(ecm_go[, c("Description","NES","pvalue","p.adjust")], row.names = FALSE)
  }
}

# 画图（TOP15）
p_go <- dotplot(gse_go, showCategory = 15,
                title = "GREM1 Knockout - GSEA GO (BP)",
                font.size = 10)
ggsave(file.path(figure_dir, "10_GREM1_GSEA_GO.png"), p_go,
       width = 10, height = 8)

write.csv(as.data.frame(gse_go), file.path(result_dir, "GREM1_GSEA_GO.csv"),
          row.names = FALSE)

# ================================================================
# 5. KEGG GSEA 富集
# ================================================================
cat("\n=== KEGG GSEA 富集 ===\n")
gse_kegg <- gseKEGG(geneList = gene_list, pvalueCutoff = 0.5,
                    verbose = FALSE)

cat("KEGG 通路数 (p<0.5):", nrow(gse_kegg), "\n")

gse_kegg_sig <- as.data.frame(gse_kegg) %>%
  filter(p.adjust < 0.25) %>%
  arrange(desc(NES))

cat("KEGG 显著通路数 (p.adj<0.25):", nrow(gse_kegg_sig), "\n")

if (nrow(gse_kegg_sig) > 0) {
  cat("\n--- TOP15 KEGG 通路 ---\n")
  print(head(gse_kegg_sig[, c("Description","NES","pvalue","p.adjust")], 15),
        row.names = FALSE)
}

if (nrow(gse_kegg) > 0) {
  p_kegg <- dotplot(gse_kegg, showCategory = 15,
                    title = "GREM1 Knockout - GSEA KEGG",
                    font.size = 10)
  ggsave(file.path(figure_dir, "11_GREM1_GSEA_KEGG.png"), p_kegg,
         width = 8, height = 6)
  write.csv(as.data.frame(gse_kegg), file.path(result_dir, "GREM1_GSEA_KEGG.csv"),
            row.names = FALSE)
}

# ================================================================
# 6. 保存完整结果
# ================================================================
write.csv(all_genes, file.path(result_dir, "GREM1_Knockout_all_genes.csv"),
          row.names = FALSE)


# 保存 knk_result
saveRDS(knk_result, file.path(result_dir, "scTenifoldKnk_GREM1.rds"))
cat("knk_result 已保存:", file.path(result_dir, "scTenifoldKnk_GREM1.rds"), "\n")

# 同时保存 GSEA 结果
saveRDS(gse_go, file.path(result_dir, "GSEA_GO_GREM1.rds"))
saveRDS(gse_kegg, file.path(result_dir, "GSEA_KEGG_GREM1.rds"))
cat("GSEA 结果已保存\n")

# 保存显著通路表格
write.csv(as.data.frame(gse_go), file.path(result_dir, "GREM1_GSEA_GO_all.csv"), row.names = FALSE)
write.csv(as.data.frame(gse_kegg), file.path(result_dir, "GREM1_GSEA_KEGG_all.csv"), row.names = FALSE)
cat("全部保存完成！\n")
