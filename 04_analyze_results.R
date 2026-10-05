# ================================================================
# Step 4: 虚拟筛选结果分析
# 功能：结合能排序、TOP20 筛选、氢键分析、富集分析、作图
# ================================================================

library(ggplot2)
library(dplyr)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)

# ===== 目录设置 =====
base_dir <- "D:/Bioinfor/colorectal/MR/virtual_screening"
result_dir <- file.path(base_dir, "04_results")
fig_dir <- file.path(result_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ===== 1. 读取对接结果 =====
docking <- read.csv(file.path(result_dir, "docking_scores.csv"))
cat("总化合物数:", nrow(docking), "\n")
cat("成功对接:", sum(docking$status == "success"), "\n")

# 过滤有结合能的
docking_valid <- docking %>%
  filter(!is.na(best_score)) %>%
  arrange(best_score)

cat("有效结合能:", nrow(docking_valid), "\n")

# ===== 2. 结合能分布 =====
p_dist <- ggplot(docking_valid, aes(x = best_score)) +
  geom_histogram(binwidth = 0.5, fill = "#4DBBD5", color = "white") +
  geom_vline(xintercept = -7.0, linetype = "dashed", color = "red", size = 0.8) +
  geom_vline(xintercept = -5.0, linetype = "dashed", color = "orange", size = 0.8) +
  annotate("text", x = -7.5, y = Inf, label = "强结合\n≤ -7.0",
           hjust = 0, vjust = 2, color = "red", size = 3.5) +
  annotate("text", x = -5.5, y = Inf, label = "中等结合\n≤ -5.0",
           hjust = 0, vjust = 2, color = "orange", size = 3.5) +
  labs(x = "Binding Affinity (kcal/mol)", y = "Count",
       title = "Distribution of Docking Scores") +
  theme_bw()

ggsave(file.path(fig_dir, "01_docking_score_distribution.png"),
       p_dist, width = 8, height = 5, dpi = 300)

# ===== 3. TOP20 化合物 =====
top20 <- head(docking_valid, 20)
write.csv(top20, file.path(result_dir, "top20_compounds.csv"), row.names = FALSE)

cat("\n=== TOP20 化合物 ===\n")
print(top20[, c("ligand", "best_score")], row.names = FALSE)

# TOP20 柱状图
p_top20 <- ggplot(top20, aes(x = reorder(ligand, -best_score), y = best_score)) +
  geom_bar(stat = "identity", fill = "#E64B35", width = 0.7) +
  coord_flip() +
  labs(x = "", y = "Binding Affinity (kcal/mol)",
       title = "TOP 20 Compounds by Binding Affinity") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 9))

ggsave(file.path(fig_dir, "02_top20_compounds.png"),
       p_top20, width = 8, height = 8, dpi = 300)

# ===== 4. 强结合化合物筛选（≤ -7.0）=====
strong_binders <- docking_valid %>%
  filter(best_score <= -7.0)

cat("\n强结合化合物 (≤ -7.0 kcal/mol):", nrow(strong_binders), "\n")
if (nrow(strong_binders) > 0) {
  write.csv(strong_binders, file.path(result_dir, "strong_binders.csv"), row.names = FALSE)
}

# ===== 5. 已知药物对比（如果有 Ginisortamab 等已知抑制剂的对接结果）=====
# 可以把已知抑制剂的 SMILES 也加入化合物库，作为阳性对照
# known_drugs <- c("Ginisortamab", "TST003")  # 这些是单抗，不适合小分子对接
# 小分子阳性对照可以用已知的 BMP 通路抑制剂

# ===== 6. TOP 化合物的靶点预测与富集分析 =====
# 如果有化合物的已知靶点数据（从 DrugBank/TCMSP 获取），可以做富集分析
# 这里假设已有靶点映射文件 compound_targets.csv

target_file <- file.path(base_dir, "02_ligands", "compound_targets.csv")
if (file.exists(target_file)) {
  targets <- read.csv(target_file)

  # TOP20 化合物的靶点
  top20_targets <- targets %>%
    filter(ligand %in% top20$ligand) %>%
    distinct(target)

  cat("\nTOP20 化合物靶点总数:", nrow(top20_targets), "\n")

  # ID 转换
  eg <- bitr(top20_targets$target, fromType = "SYMBOL",
             toType = "ENTREZID", OrgDb = org.Hs.eg.db)

  # GO 富集
  go_res <- enrichGO(gene = eg$ENTREZID, OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID", ont = "BP",
                     pvalueCutoff = 0.5)

  if (nrow(go_res) > 0) {
    p_go <- dotplot(go_res, showCategory = 15,
                    title = "TOP20 Compounds - Target GO Enrichment")
    ggsave(file.path(fig_dir, "03_target_GO_enrichment.png"),
           p_go, width = 10, height = 7, dpi = 300)
    write.csv(as.data.frame(go_res), file.path(result_dir, "target_GO.csv"), row.names = FALSE)
  }

  # KEGG 富集
  kegg_res <- enrichKEGG(gene = eg$ENTREZID, pvalueCutoff = 0.5)
  if (nrow(kegg_res) > 0) {
    p_kegg <- dotplot(kegg_res, showCategory = 15,
                      title = "TOP20 Compounds - Target KEGG Enrichment")
    ggsave(file.path(fig_dir, "04_target_KEGG_enrichment.png"),
           p_kegg, width = 8, height = 6, dpi = 300)
    write.csv(as.data.frame(kegg_res), file.path(result_dir, "target_KEGG.csv"), row.names = FALSE)
  }
}

# ===== 7. 化合物-靶点-通路网络（Cytoscape 输入）=====
# 生成 Cytoscape 可用的边表和节点表
if (file.exists(target_file)) {
  network_edges <- targets %>%
    filter(ligand %in% top20$ligand) %>%
    select(source = ligand, target = target) %>%
    mutate(interaction = "targets")

  write.csv(network_edges, file.path(result_dir, "network_edges.csv"), row.names = FALSE)
  cat("\nCytoscape 网络边表已保存\n")
}

# ===== 8. 汇总报告 =====
sink(file.path(result_dir, "docking_summary.txt"))
cat("========================================\n")
cat("GREM1 虚拟筛选结果汇总\n")
cat("========================================\n\n")
cat("受体: PDB 5AEJ (GREM1)\n")
cat("化合物库: DrugBank + TCMSP\n")
cat("对接软件: AutoDock Vina\n")
cat("对接盒子: 20x20x20 A\n\n")
cat("总化合物数:", nrow(docking), "\n")
cat("成功对接:", sum(docking$status == "success"), "\n")
cat("强结合 (≤ -7.0):", nrow(strong_binders), "\n")
cat("中等结合 (-7.0 ~ -5.0):", sum(docking_valid$best_score > -7.0 & docking_valid$best_score <= -5.0), "\n\n")
cat("=== TOP10 化合物 ===\n")
print(head(docking_valid[, c("ligand", "best_score")], 10), row.names = FALSE)
cat("\n========================================\n")
sink()

cat("\n========================================\n")
cat("分析完成！结果保存在:", result_dir, "\n")
cat("  - docking_scores.csv (全部结合能)\n")
cat("  - top20_compounds.csv (TOP20)\n")
cat("  - strong_binders.csv (强结合)\n")
cat("  - figures/ (图表)\n")
cat("  - docking_summary.txt (汇总报告)\n")
cat("========================================\n")
