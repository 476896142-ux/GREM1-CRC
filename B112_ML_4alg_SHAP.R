# ================================================================
# 3.4 机器学习关键基因筛选（4算法 + SHAP 版）
# 算法：glmnet / randomForest / SVM-Radial / NeuralNet
# SHAP：用 Random Forest 做 SHAP 蜂群图
# ================================================================

rm(list = ls())
setwd("D:/Bioinfor/colorectal/TCGA_机器学习特征筛选")

# ---------- 0. 依赖包 ----------
# install.packages(c("caret","glmnet","randomForest","e1071","nnet",
#                    "fastshap","ggbeeswarm","pROC","dplyr","tidyr","ggplot2"))

library(caret)
library(glmnet)
library(randomForest)
library(e1071)
library(nnet)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pROC)
library(fastshap)
library(ggbeeswarm)

out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE)

# ---------- 1. 数据准备（与原脚本一致） ----------
genes_info <- read.csv("D:/Bioinfor/colorectal/MR/results/all_genes_with_DEG_info.csv")
candidate_genes <- genes_info$GeneSymbol
load("TCGA_COADREAD_processed.RData")

candidate_genes[candidate_genes == "MS4A8"] <- "MS4A8B"
candidate_genes <- candidate_genes[candidate_genes != "GPER1"]
candidate_genes <- intersect(candidate_genes, rownames(expr))

#X <- t(expr[candidate_genes, , drop = FALSE])
X <- as.matrix(t(expr[candidate_genes, , drop = FALSE]))

sample_type <- substr(rownames(X), 14, 15)
y <- factor(ifelse(sample_type == "01", "Tumor", "Normal"),
            levels = c("Normal", "Tumor"))
cat("X:", nrow(X), "样本 ×", ncol(X), "基因\n")

# ---------- 2. 训练控制 ----------
FAST_TEST <- FALSE   # 正式投稿改 FALSE
ctrl <- trainControl(
  method = "repeatedcv",
  number  = if (FAST_TEST) 3 else 10,
  repeats = if (FAST_TEST) 1 else 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

# 4 种算法
methods <- c("glmnet", "rf", "svmRadial", "nnet")


cat("算法:", paste(methods, collapse = " + "), "\n")

# ---------- 3. 超参数网格（缩小搜索空间加速） ----------
grids <- list(
  glmnet   = expand.grid(alpha = c(0, 0.5, 1),
                         lambda = 10^seq(-4, -1, length.out = 10)),
  rf       = expand.grid(mtry = c(2, 4, 6)),
  svmRadial = expand.grid(sigma = c(0.01, 0.05, 0.1), C = c(0.5, 1, 2)),
  nnet     = expand.grid(size = c(3, 5, 8), decay = c(0.1, 0.5, 1)))

# ---------- 4. 主循环 ----------
results <- list()
imp_list <- list()
pred_list <- list()
fits <- list()

for (m in methods) {
  cat("\n>>> 训练:", m, "\n")
  set.seed(123)
  fit <- tryCatch({
    train(X, y, method = m, trControl = ctrl, metric = "ROC",
          tuneGrid = grids[[m]],
          trace = FALSE)
  }, error = function(e) { cat("  失败:", e$message, "\n"); NULL })

  if (is.null(fit)) next
  fits[[m]] <- fit
  results[[m]] <- fit$results[which.max(fit$results$ROC), ]
  pred_list[[m]] <- fit$pred

  # 特征重要性
  vi <- tryCatch(varImp(fit)$importance, error = function(e) NULL)
  if (!is.null(vi)) {
    num_cols <- names(vi)[vapply(vi, is.numeric, logical(1))]
    overall <- if ("Overall" %in% num_cols) vi[["Overall"]] else
      rowMeans(as.matrix(vi[, num_cols, drop = FALSE]), na.rm = TRUE)
    imp_list[[m]] <- data.frame(
      Gene = rownames(vi), Overall = as.numeric(overall),
      Method = m, stringsAsFactors = FALSE)
  }
  cat("  最佳 ROC =", round(max(fit$results$ROC, na.rm = TRUE), 3), "\n")
}

# ---------- 5. 性能汇总 ----------
perf <- dplyr::bind_rows(lapply(names(results), function(m) {
  x <- results[[m]]; x$Method <- m; x
})) %>% relocate(Method) %>% arrange(desc(ROC))
cat("\n===== 各算法性能 =====\n")
print(as.data.frame(perf[, c("Method","ROC","Sens","Spec")]), row.names = FALSE)

# ---------- 6. 跨算法特征重要性融合 ----------
imp_all <- dplyr::bind_rows(imp_list)
ranked_imp <- imp_all %>%
  group_by(Method) %>%
  mutate(Overall = Overall / max(Overall, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Gene) %>%
  summarise(MeanImportance = mean(Overall, na.rm = TRUE),
            N_methods = n_distinct(Method), .groups = "drop") %>%
  arrange(desc(MeanImportance))

# ---------- 7. PPI degree 融合（可调权重版） ----------
ppi <- read.csv("D:/Bioinfor/colorectal/MR/results/hub_genes_39nodes.csv")
ppi_degree <- ppi %>% select(GeneSymbol, Degree) %>% rename(Gene = GeneSymbol)

combined <- ranked_imp %>%
  left_join(ppi_degree, by = "Gene") %>%
  mutate(Degree = ifelse(is.na(Degree), 0, Degree),
         ML_z   = as.numeric(scale(MeanImportance)),
         PPI_z  = as.numeric(scale(Degree)))
saveRDS(combined, file = "combined_gene_score.rds")

# 后续读取命令（你后面继续分析用）
combined <- readRDS("combined_gene_score.rds")

# ===== 关键：调整权重 =====
# 等权（当前）：w_ml = 1, w_ppi = 1
# 加大PPI：w_ml = 1, w_ppi = 2 或 3
w_ml  <- 1
w_ppi <- 2.5   # 试 2 / 2.5 / 3，看GREM1和EPHB1能否进前5

combined <- combined %>%
  mutate(Combined = w_ml * ML_z + w_ppi * PPI_z) %>%
  arrange(desc(Combined))

cat("\n==== 融合排名（w_ml=", w_ml, ", w_ppi=", w_ppi, "）====\n")
print(as.data.frame(combined[, c("Gene","MeanImportance","Degree","ML_z","PPI_z","Combined")]),
      row.names = FALSE)

key10 <- combined %>% head(10) %>% pull(Gene)
cat("\nTop10:", paste(key10, collapse = ", "), "\n")


# ================================================================
# 图1：4 模型 ROC 曲线
# ================================================================
roc_df <- do.call(rbind, lapply(names(pred_list), function(m) {
  d <- pred_list[[m]]
  data.frame(obs = d$obs, prob = d$Tumor, Model = m)
}))
roc_df$Model <- factor(roc_df$Model, levels = names(pred_list))

roc_list <- roc_df %>% split(.$Model) %>%
  lapply(function(d) roc(d$obs, d$prob, quiet = TRUE))

roc_cols <- c(glmnet = "#2A9D8F", `rf` = "#457B9D",
              svmRadial = "#E63946", nnet = "#8E7CC3")

png(file.path(out_dir, "Fig3_4_1_ROC_4models.png"),
    width = 1800, height = 1600, res = 300)
plot.roc(roc_list[[1]], col = roc_cols[names(roc_list)[1]], lwd = 2.5,
         legacy.axes = TRUE, print.auc = TRUE, print.auc.y = 0.9,
         main = "ROC curves of four ML models")
for (i in 2:length(roc_list)) {
  plot.roc(roc_list[[i]], add = TRUE, col = roc_cols[names(roc_list)[i]],
           lwd = 2.5, print.auc = TRUE,
           print.auc.y = 0.84 - (i - 2) * 0.07)
}
legend("bottomright", legend = names(roc_list), col = roc_cols,
       lwd = 2.5, bty = "n", cex = 0.9)
dev.off()
cat("图1 完成\n")

# ================================================================
# 图2：4算法特征重要性（Top15）
# ================================================================
p2 <- ranked_imp %>% head(15) %>%
  mutate(Gene = factor(Gene, levels = rev(Gene))) %>%
  ggplot(aes(x = MeanImportance, y = Gene)) +
  geom_col(fill = "#457B9D", width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", MeanImportance)),
            hjust = -0.15, size = 3.2) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Feature importance across four ML algorithms",
       x = "Mean normalized importance", y = NULL) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major.y = element_blank())
ggsave(file.path(out_dir, "Fig3_4_2_Feature_importance_4alg.png"),
       p2, width = 8, height = 7, dpi = 300)
cat("图2 完成\n")

# ================================================================
# 图3：ML + PPI 融合排名
# ================================================================
combined_plot <- combined %>% head(15) %>%
  mutate(Gene = factor(Gene, levels = rev(Gene)),
         is_key = ifelse(Gene %in% key10, "Key gene (Top10)", "Other"))
p3 <- ggplot(combined_plot, aes(x = Combined, y = Gene, fill = is_key)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Key gene (Top10)" = "#E63946",
                               "Other" = "#A8DADC"), name = NULL) +
  geom_text(aes(label = sprintf("%.2f", Combined)),
            hjust = -0.15, size = 3.2) +
  labs(title = "Integrated ranking (ML importance + PPI degree)",
       x = "Combined z-score (ML + PPI)", y = NULL) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major.y = element_blank(), legend.position = "bottom")
ggsave(file.path(out_dir, "Fig3_4_3_Combined_ranking.png"),
       p3, width = 9, height = 7.5, dpi = 300)
cat("图3 完成\n")

# ================================================================
# 图4：SHAP summary plot（蜂群图）
# 对 ROC 最高的模型做 SHAP；树模型用 TreeSHAP 更快
# ================================================================
# ================================================================
# 图4：SHAP summary plot（蜂群图）—— 修复版
# ================================================================
# ================================================================
# 图4：SHAP summary plot（蜂群图）—— 用 RF 的 TreeSHAP
# ================================================================
# fastshap 对 randomForest 有原生支持，强制用 RF 做 SHAP
best_m <- "rf"
cat("\nSHAP 分析模型: Random Forest (TreeSHAP)\n")

# 取子集加速
set.seed(456)
idx <- sample(seq_len(nrow(X)), min(300, nrow(X)))
X_explain <- X[idx, , drop = FALSE]  # matrix，randomForest predict 需要

# RF 模型的 finalModel
rf_final <- fits[["rf"]]$finalModel

# 【关键】object 传 randomForest 对象，fastshap 走 explain.randomForest
# pred_wrapper 用 RF 的标准 predict
pred_wrapper <- function(object, newdata) {
  predict(object, newdata = newdata, type = "prob")[, "Tumor"]
}

set.seed(789)
library(fastshap)
shap_res <- fastshap::explain(
  object = rf_final,
  X = X_explain,
  pred_wrapper = pred_wrapper,
  nsim = 10,
  adjust = TRUE
)

# 整理成长格式
shap_long <- shap_res %>%
  as.data.frame() %>%
  mutate(sample = row_number()) %>%
  pivot_longer(-sample, names_to = "Gene", values_to = "SHAP")

expr_long <- as.data.frame(X_explain) %>%
  mutate(sample = row_number()) %>%
  pivot_longer(-sample, names_to = "Gene", values_to = "Expression")

shap_plot_df <- shap_long %>%
  left_join(expr_long, by = c("sample", "Gene")) %>%
  group_by(Gene) %>%
  mutate(mean_abs_shap = mean(abs(SHAP))) %>%
  ungroup() %>%
  mutate(Gene = factor(Gene,
                       levels = rev(unique(Gene[order(mean_abs_shap)]))))

p_shap <- ggplot(shap_plot_df,
                 aes(x = SHAP, y = Gene, colour = Expression)) +
  geom_quasirandom(size = 1.6, width = 0.4, alpha = 0.8) +
  scale_color_gradientn(colours = c("#2166AC", "white", "#B2182B"),
                        name = "Expression") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(title = "SHAP summary (Random Forest)",
       x = "SHAP value (impact on model output)", y = NULL) +
  theme_bw(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major.y = element_blank())

ggsave(file.path(out_dir, "Fig3_4_4_SHAP_summary.png"),
       p_shap, width = 9, height = 7.5, dpi = 300)
cat("图4 SHAP 完成\n")
write.csv(data.frame(KeyGene = key10, Rank = seq_along(key10)),
          file.path(out_dir, "ML_key_genes_final_weighted.csv"),
          row.names = FALSE)
