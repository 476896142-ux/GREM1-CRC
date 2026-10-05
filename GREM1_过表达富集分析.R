# ================================================================
# GREM1 过表达响应基因 GO/KEGG 富集分析
# 输入：GREM1_Overexpression_TOP20.csv
# 输出：GO/KEGG 富集表 + 气泡图
# ================================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# ===== 目录设置 =====
result_dir <- "D:/Bioinfor/colorectal/MR/results/single_cell"
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

# ================================================================
# 1. 读取过表达 TOP20 响应基因
# ================================================================
top20_oe <- read.csv(file.path(result_dir, "GREM1_Overexpression_TOP20.csv"))
cat("过表达 TOP20 基因:", paste(top20_oe$gene, collapse = ", "), "\n")

# ================================================================
# 2. ID 转换（SYMBOL → ENTREZID）
# ================================================================
eg_oe <- bitr(top20_oe$gene, fromType = "SYMBOL",
              toType = "ENTREZID", OrgDb = org.Hs.eg.db)
cat("ID 映射成功:", nrow(eg_oe), "/ 20\n")

# ================================================================
# 3. GO 富集分析（BP/CC/MF 全部）
# ================================================================
go_oe <- enrichGO(gene = eg_oe$ENTREZID, OrgDb = org.Hs.eg.db,
                  keyType = "ENTREZID", ont = "ALL",
                  pvalueCutoff = 1, qvalueCutoff = 1)

go_oe_df <- as.data.frame(go_oe) %>%
  arrange(pvalue) %>%
  head(15)

cat("\n=== 过表达 TOP20 GO 富集 TOP15 ===\n")
print(go_oe_df[, c("Description","ONTOLOGY","pvalue","p.adjust")], row.names = FALSE)

# ================================================================
# 4. KEGG 富集分析
# ================================================================
kegg_oe <- enrichKEGG(gene = eg_oe$ENTREZID,
                      pvalueCutoff = 1, qvalueCutoff = 1)

kegg_oe_df <- as.data.frame(kegg_oe) %>%
  arrange(pvalue) %>%
  head(15)

cat("\n=== 过表达 TOP20 KEGG 富集 TOP15 ===\n")
print(kegg_oe_df[, c("Description","pvalue","p.adjust")], row.names = FALSE)

# ================================================================
# 5. 保存富集结果
# ================================================================
write.csv(as.data.frame(go_oe),
          file.path(result_dir, "GREM1_OE_GO_ORA.csv"),
          row.names = FALSE)
write.csv(as.data.frame(kegg_oe),
          file.path(result_dir, "GREM1_OE_KEGG_ORA.csv"),
          row.names = FALSE)

# ================================================================
# 6. 绘制气泡图
# ================================================================
if (nrow(go_oe_df) > 0) {
  p_go <- dotplot(go_oe, showCategory = 10) +
    facet_grid(ONTOLOGY ~ ., scales = "free_y") +
    ggtitle("GREM1 Overexpression TOP20 - GO Enrichment") +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file.path(result_dir, "12_GREM1_OE_GO.png"),
         p_go, width = 10, height = 8, dpi = 300)
  cat("\nGO 气泡图已保存: 12_GREM1_OE_GO.png\n")
}

if (nrow(kegg_oe_df) > 0) {
  p_kegg <- dotplot(kegg_oe, showCategory = 10) +
    ggtitle("GREM1 Overexpression TOP20 - KEGG") +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file.path(result_dir, "13_GREM1_OE_KEGG.png"),
         p_kegg, width = 8, height = 6, dpi = 300)
  cat("KEGG 气泡图已保存: 13_GREM1_OE_KEGG.png\n")
}

# ================================================================
# 7. 保存 R 对象（方便后续读取）
# ================================================================
saveRDS(go_oe, file.path(result_dir, "GREM1_OE_GO_ORA.rds"))
saveRDS(kegg_oe, file.path(result_dir, "GREM1_OE_KEGG_ORA.rds"))

cat("\n========================================\n")
cat("过表达富集分析完成！\n")
cat("结果保存在:", result_dir, "\n")
cat("  - GREM1_OE_GO_ORA.csv\n")
cat("  - GREM1_OE_KEGG_ORA.csv\n")
cat("  - 12_GREM1_OE_GO.png\n")
cat("  - 13_GREM1_OE_KEGG.png\n")
cat("========================================\n")
