# ===== B2: 一致性聚类分亚型 =====
BiocManager::install("ConsensusClusterPlus")
library(ConsensusClusterPlus)
key_genes_final_df<-read.csv("TCGA_机器学习特征筛选/ML_key_genes_final_weighted.csv")
key_genes_final <- key_genes_final_df$KeyGene
print(key_genes_final)
# 仅肿瘤样本的 key genes 表达
tumor_idx <- rownames(X)[y == "Tumor"]
class(X)
str(X, max.level = 1)

expr_key <- t(X[tumor_idx, key_genes_final, drop = FALSE])   # 基因 × 肿瘤样本


# 一致性聚类（k = 2~6）
cc <- ConsensusClusterPlus(
  as.matrix(expr_key),
  maxK = 6, reps = 1000, pItem = 0.8, pFeature = 1,
  clusterAlg = "pam", distance = "pearson",
  seed = 123, plot = "png",
  title = "MR/results/consensus_cluster"
)

# 看 CDF 曲线选最优 k（通常 k=2 或 3）
# 提取亚型标签（选 k=4）
subtype <- as.data.frame(cc[[4]]$consensusClass)
colnames(subtype) <- "Subtype"
subtype$sample <- rownames(subtype)
table(subtype$Subtype)

write.csv(subtype, "MR/results/consensus_subtype.csv", row.names = FALSE)


# 计算每个 K 的 PAC 值
calc_pac <- function(consensus_matrix) {
  cm <- consensus_matrix[upper.tri(consensus_matrix)]
  sum(cm > 0.1 & cm < 0.9) / length(cm)
}

pac_values <- sapply(2:6, function(k) calc_pac(cc[[k]]$consensusMatrix))
names(pac_values) <- paste0("k=", 2:6)
print(pac_values)
cat("PAC 最小的 K:", names(which.min(pac_values)), "\n")
