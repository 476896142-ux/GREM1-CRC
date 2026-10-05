# 4. 9 个 key genes 细胞定位
# ================================================================
key_genes <- c("GREM1","CSE1L","KLB","ANGPTL1","UGDH",
               "EPB41L3","ZNF239","NTRK3","ITIH5")
genes_in_data <- key_genes[key_genes %in% rownames(sc)]
cat("在数据中找到的基因:", paste(genes_in_data, collapse = ", "), "\n")

# 每个基因的 UMAP 表达图
for (g in genes_in_data) {
  p <- FeaturePlot(sc, features = g, reduction = "umap") +
    ggtitle(g) + theme(plot.title = element_text(hjust = 0.5))
  ggsave(file.path(figure_dir, paste0("03_expr_", g, ".png")), p, width = 7, height = 6)
}

# 多基因联合 UMAP 图（重点基因）
top_genes <- intersect(c("GREM1","CSE1L","NTRK3","LAMC2"), rownames(sc))
if (length(top_genes) >= 2) {
  p_multi <- FeaturePlot(sc, features = top_genes, reduction = "umap", ncol = 2)
  ggsave(file.path(figure_dir, "04_expr_key_genes_combined.png"),
         p_multi, width = 12, height = 10)
}

# 各细胞类型平均表达热图
avg_expr <- AverageExpression(sc, features = genes_in_data, group.by = "celltype")[[1]]
avg_expr_scaled <- t(scale(t(avg_expr)))
library(pheatmap)
pheatmap(avg_expr_scaled,
         main = "Key Genes Expression by Cell Type (z-score)",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         filename = file.path(figure_dir, "05_keygenes_heatmap.png"),
         width = 10, height = 6)

# GREM1 各细胞类型表达详情
grem1_detail <- FetchData(sc, vars = c("GREM1","celltype")) %>%
  group_by(celltype) %>%
  summarise(mean_expr = mean(GREM1),
            pct_expressing = round(mean(GREM1 > 0) * 100, 1),
            n_cells = n(), .groups = "drop") %>%
  arrange(desc(mean_expr))
cat("\n=== GREM1 细胞类型表达 ===\n")
print(as.data.frame(grem1_detail), row.names = FALSE)
write.csv(grem1_detail, file.path(result_dir, "GREM1_expression_by_celltype.csv"),
          row.names = FALSE)
# 5. CellChat 细胞通讯
# ================================================================
library(CellChat)
# ===== 检查并清理 celltype =====
table(sc$celltype, useNA = "ifany")   # 看有没有 NA

# 去掉 celltype 为 NA 的细胞
sc <- subset(sc, cells = colnames(sc)[!is.na(sc$celltype) & sc$celltype != ""])
sc$celltype <- droplevels(as.factor(sc$celltype))
cat("清理后细胞数:", ncol(sc), "\n")
cat("细胞类型:", paste(levels(sc$celltype), collapse = ", "), "\n")

# ===== 重新提取并创建 CellChat =====
expr_matrix <- GetAssayData(sc, layer = "data")
cell_meta   <- sc@meta.data

cellchat <- createCellChat(object = expr_matrix, meta = cell_meta, group.by = "celltype")

# 确认 idents 正确
print(unique(cellchat@idents))

# ===== 后续正常跑 =====
cellchat@DB <- CellChatDB.human
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

cat("CellChat 完成！\n")
# ===== 保存 + 画图 =====
saveRDS(sc, file.path(result_dir, "seurat_merged.rds"))
saveRDS(cellchat, file.path(result_dir, "cellchat_object.rds"))

png(file.path(figure_dir, "06_cellchat_overall.png"), width = 1000, height = 800, res = 150)
netVisual_circle(cellchat@net$count,
                 vertex.weight = rowSums(cellchat@net$count),
                 weight.scale = TRUE, label.edge = FALSE,
                 vertex.label.cex = 0.8)
dev.off()
# BMP 通路
bmp_pathways <- grep("BMP", cellchat@netP$pathways, value = TRUE)
cat("BMP 通路:", paste(bmp_pathways, collapse = ", "), "\n")
for (pw in bmp_pathways) {
  png(file.path(result_dir, paste0("07_cellchat_", pw, ".png")),
      width = 900, height = 700, res = 150)
  netVisual_aggregate(cellchat, signaling = pw, layout = "circle")
  dev.off()
}

# 配体-受体对
lr_pairs <- subsetCommunication(cellchat)
write.csv(lr_pairs, file.path(result_dir, "cellchat_lr_pairs.csv"), row.names = FALSE)

# GREM1 相关
grem1_lr <- lr_pairs[grep("GREM1", lr_pairs$ligand) | grep("GREM1", lr_pairs$receptor), ]
if (nrow(grem1_lr) > 0) {
  cat("\n=== GREM1 相关配体-受体对 ===\n")
  print(as.data.frame(grem1_lr[, c("source","target","ligand","receptor","prob")]), row.names = FALSE)
  write.csv(grem1_lr, file.path(result_dir, "GREM1_lr_pairs.csv"), row.names = FALSE)
} else {
  cat("\nGREM1 不在 CellChat 配体-受体数据库中（它是分泌型 BMP 拮抗剂，可能以非经典方式作用）\n")
}
