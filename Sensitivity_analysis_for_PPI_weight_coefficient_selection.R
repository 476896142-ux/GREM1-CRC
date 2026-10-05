library(dplyr)
library(tibble)
library(purrr)
library(stringr)

##==================== 1. 基础数据预处理（沿用你前面的代码）====================
# 你的原始代码
combined_base <- ranked_imp %>%
  left_join(ppi_degree, by = "Gene") %>%
  mutate(
    Degree = ifelse(is.na(Degree), 0, Degree),
    ML_z   = as.numeric(scale(MeanImportance)),
    PPI_z  = as.numeric(scale(Degree))
  )

# 待测试权重梯度
weight_vec <- c(1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0)

##==================== 2. 循环批量计算每个权重、提取top10 ====================
# 构建列表存储结果
sensitivity_res <- list()
top10_gene_list <- list()

for(w in weight_vec){
  df <- combined_base %>%
    mutate(Combined_score = ML_z + w * PPI_z) %>%
    arrange(desc(Combined_score))
  
  # 保存全部打分结果
  sensitivity_res[[as.character(w)]] <- df
  
  # 提取top10基因
  top10 <- df %>% slice_head(n = 10)
  top10_gene_list[[as.character(w)]] <- top10$Gene
  
  # 输出汇总指标：平均PPI degree、平均ML_z
  metrics <- top10 %>%
    summarise(
      weight = w,
      mean_degree = mean(Degree),
      mean_ML_z = mean(ML_z),
      mean_PPI_z = mean(PPI_z),
      mean_Combined = mean(Combined_score)
    )
  
  if(!exists("metric_df")){
    metric_df <- metrics
  }else{
    metric_df <- bind_rows(metric_df, metrics)
  }
}

# 保存指标表（用于折线图）
write.csv(metric_df, "sensitivity_metrics.csv", row.names = FALSE)

# 保存每个权重下top10基因
saveRDS(top10_gene_list, "sensitivity_top10_genelist.rds")

##==================== 3. 计算Jaccard相似度矩阵 ====================
# Jaccard = 交集基因数 / 并集基因数
jaccard_mat <- matrix(0, nrow = length(weight_vec), ncol = length(weight_vec),
                      dimnames = list(as.character(weight_vec), as.character(weight_vec)))

for(i in seq_along(weight_vec)){
  for(j in seq_along(weight_vec)){
    g1 <- top10_gene_list[[i]]
    g2 <- top10_gene_list[[j]]
    inter <- length(intersect(g1,g2))
    union <- length(union(g1,g2))
    jaccard_mat[i,j] <- inter / union
  }
}

# 保存Jaccard矩阵
write.csv(jaccard_mat, "sensitivity_jaccard_matrix.csv", row.names = TRUE)
saveRDS(jaccard_mat, "sensitivity_jaccard_matrix.rds")

##==================== 4. 绘图：A Jaccard热图 + B 双指标折线图 ====================
library(ggplot2)
library(patchwork)
library(reshape2)

# 图A：Jaccard热图
jacc_long <- melt(jaccard_mat)
p_heat <- ggplot(jacc_long, aes(Var1, Var2, fill=value)) +
  geom_tile() +
  scale_fill_viridis_c(limits = c(0,1)) +
  labs(x="PPI weight", y="PPI weight", fill="Jaccard similarity") +
  theme_bw() +
  theme(aspect.ratio = 1)

# 图B：双轴折线图：mean_degree 和 mean_ML_z
p_line <- ggplot(metric_df, aes(x = weight)) +
  geom_line(aes(y = mean_degree, color="Mean PPI Degree"), linewidth=1) +
  geom_point(aes(y = mean_degree, color="Mean PPI Degree"), size=2) +
  geom_line(aes(y = mean_ML_z*25, color="Mean ML_z"), linewidth=1) +
  geom_point(aes(y = mean_ML_z*25, color="Mean ML_z"), size=2) +
  geom_vline(xintercept = 2.5, linetype="dashed", color="red", linewidth=0.8)+
  scale_y_continuous(name="Mean PPI Degree",
                     sec.axis = sec_axis(~./25, name="Mean ML_z")) +
  scale_color_manual(values = c("#2E86AB","#A23B72"))+
  labs(x="PPI weight coefficient", color="Metrics") +
  theme_bw()

# 拼接图A+B
p_combined <- p_heat + p_line + plot_annotation(tag_levels = "A")
ggsave("Sensitivity_figure.pdf", p_combined, width = 10, height = 4.5)
ggsave("Sensitivity_figure.png", p_combined, width = 10, height = 4.5, dpi=300)

library(clusterProfiler)
library(org.Hs.eg.db)

enrich_count <- c()
for(w in names(top10_gene_list)){
  genes <- top10_gene_list[[w]]
  ego <- enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                  ont = "BP", pAdjustMethod = "fdr", qvalueCutoff = 0.05)
  cnt <- nrow(ego@result %>% filter(p.adjust < 0.05))
  enrich_count <- c(enrich_count, cnt)
}
metric_df$enrich_term_count <- enrich_count

# 富集条目柱状图
p_bar <- ggplot(metric_df, aes(x=weight, y=enrich_term_count))+
  geom_col(fill="#347eb8") +
  geom_vline(xintercept = 2.5, linetype="dashed",color="red")+
  labs(x="PPI weight coefficient",y="Number of enriched terms (FDR<0.05)")+
  theme_bw()

# 重新拼接成4-panel图
p_heat + p_line + p_bar

library(patchwork)
library(ggplot2)

# 拼接三张图，按 A/B/C 自动编号
p_4panel <- p_heat + p_line + p_bar +
  plot_annotation(
    tag_levels = "A",
    tag_suffix = ""
  ) +
  plot_layout(ncol = 3)  # 一行三列，根据你的版式也可以改成 nrow = 2

# 方法1：保存 PDF（论文推荐，矢量图）
ggsave(
  filename = "Sensitivity_4panel.pdf",
  plot = p_4panel,
  width = 14,    # 根据需要调整总宽度
  height = 4.5,  # 根据需要调整总高度
  useDingbats = FALSE
)

# 方法2：保存 PNG（高清预览）
ggsave(
  filename = "Sensitivity_4panel.png",
  plot = p_4panel,
  width = 14,
  height = 4.5,
  dpi = 300
)

# 方法3：如果想更精细控制尺寸，也可以用以下参数
# ggsave("Sensitivity_4panel.pdf", plot = p_4panel,
#        width = 14, height = 4.5, device = cairo_pdf)

