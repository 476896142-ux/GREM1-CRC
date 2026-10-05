# ================================================================

rm(list = ls())
setwd("D:/Bioinfor/colorectal")

# ---------- 1. 加载包 ----------
library(caret)
library(glmnet)
library(randomForest)
library(e1071)
library(dplyr)
library(tibble)

# ---------- 2. 数据准备 ----------
# 2.1 基因列表
genes_info <- read.csv("MR/results/all_genes_with_DEG_info.csv")
candidate_genes <- genes_info$GeneSymbol

# 2.2 TCGA 表达矩阵 + 临床
load("TCGA_机器学习特征筛选/TCGA_COADREAD_processed.RData")   # 含 expr, surv, clin, tumor_samples...

# 2.3 基因名修正
candidate_genes[candidate_genes == "MS4A8"] <- "MS4A8B"   # 表达矩阵里叫 MS4A8B
candidate_genes <- candidate_genes[candidate_genes != "GPER1"]  # 表达矩阵里没有
candidate_genes <- intersect(candidate_genes, rownames(expr))
cat("最终基因数:", length(candidate_genes), "\n")

# 2.4 构建 X（样本 × 基因）
X <- t(expr[candidate_genes, , drop = FALSE])
cat("X:", nrow(X), "样本 ×", ncol(X), "基因\n")

# 2.5 构建 y（肿瘤 vs 正常二分类）
sample_type <- substr(rownames(X), 14, 15)
y <- factor(ifelse(sample_type == "01", "Tumor", "Normal"),
            levels = c("Normal", "Tumor"))
table(y)
stopifnot(nrow(X) == length(y), !any(is.na(y)), !any(is.na(X)))

# ---------- 3. 训练控制 ----------
# 先 FAST_TEST=TRUE 跑通（3折×1次），确认后改 FALSE（10折×10次）
FAST_TEST <- TRUE

ctrl <- trainControl(
  method = "repeatedcv",
  number  = if (FAST_TEST) 3 else 10,
  repeats = if (FAST_TEST) 1 else 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)
cat("重抽样:", if (FAST_TEST) "3折×1次" else "10折×10次", "\n")

# ---------- 4. 3 种算法 ----------
methods <- c("glmnet", "rf", "svmRadial")
cat("算法:", paste(methods, collapse = " + "), "\n")

# ---------- 5. 主循环（每个算法独立 tryCatch） ----------
results <- list()
imp_list <- list()

for (m in methods) {
  cat("\n正在跑:", m, "\n")
  set.seed(123)
  fit <- tryCatch({
    train(X, y, method = m, trControl = ctrl, metric = "ROC")
  }, error = function(e) { cat("  失败:", e$message, "\n"); NULL })
  
  if (!is.null(fit)) {
    results[[m]] <- fit$results[which.max(fit$results$ROC), ]
    
    # 特征重要性：统一成 Gene + Overall 两列
    vi <- tryCatch(varImp(fit)$importance, error = function(e) NULL)
    if (!is.null(vi)) {
      if (!is.data.frame(vi)) vi <- as.data.frame(vi)
      num_cols <- names(vi)[vapply(vi, is.numeric, logical(1))]
      if (length(num_cols) > 0) {
        overall <- if ("Overall" %in% num_cols) {
          vi[["Overall"]]                            # rf/glmnet 有 Overall
        } else {
          rowMeans(as.matrix(vi[, num_cols, drop = FALSE]), na.rm = TRUE)  # svmRadial 取各类别均值
        }
        imp_list[[m]] <- data.frame(
          Gene = rownames(vi),
          Overall = as.numeric(overall),
          Method = m,
          stringsAsFactors = FALSE
        )
      }
    }
    cat("  最佳 ROC =", round(max(fit$results$ROC, na.rm = TRUE), 3), "\n")
  }
}


# ---------- 6. 模型性能汇总 ----------
perf <- dplyr::bind_rows(lapply(names(results), function(m) {
  x <- results[[m]]; x$Method <- m; x
})) %>% relocate(Method) %>% arrange(desc(ROC))

cat("\n===== 各算法性能 =====\n")
print(as.data.frame(perf[, c("Method", "ROC", "Sens", "Spec")]), row.names = FALSE)

# ---------- 7. 特征重要性综合排名 ----------
imp_all <- dplyr::bind_rows(imp_list)

# 按算法归一化到 0-1（不同算法重要性尺度不同），再取平均
ranked_imp <- imp_all %>%
  group_by(Method) %>%
  mutate(Overall = Overall / max(Overall, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(Gene) %>%
  summarise(MeanImportance = mean(Overall, na.rm = TRUE),
            N_methods = n_distinct(Method), .groups = "drop") %>%
  arrange(desc(MeanImportance))

cat("\n===== 特征重要性排名（全部）=====\n")
print(as.data.frame(ranked_imp), row.names = FALSE)

# ---------- 8. key genes：取 3 种算法共同认可的（交集） ----------
# 按 MeanImportance 排名取前 8 作为最终 key genes
key_genes_final <- ranked_imp %>%
  head(8) %>%
  pull(Gene)

cat("===== 最终 Key Genes (Top 8) =====\n")
print(as.data.frame(head(ranked_imp, 8)), row.names = FALSE)
cat("\nKey genes:", paste(key_genes_final, collapse = ", "), "\n")

######
# PPI degree（从之前的 hub_genes_39nodes.csv 读取）
ppi <- read.csv("MR/results/hub_genes_39nodes.csv")
ppi_degree <- ppi %>% select(GeneSymbol, Degree) %>%
  rename(Gene = GeneSymbol)

# 综合打分 = ML重要性 + PPI degree（各归一化后等权相加）
combined <- ranked_imp %>%
  left_join(ppi_degree, by = "Gene") %>%
  mutate(Degree = ifelse(is.na(Degree), 0, Degree),
         ML_z  = as.numeric(scale(MeanImportance)),
         PPI_z = as.numeric(scale(Degree)),
         Combined = ML_z + PPI_z) %>%
  arrange(desc(Combined))

print(as.data.frame(combined[, c("Gene","MeanImportance","Degree","Combined")]), row.names = FALSE)

# 取综合排名前 10
key_genes_final <- combined %>% head(10) %>% pull(Gene)
cat("综合打分 Top 10:", paste(key_genes_final, collapse = ", "), "\n")
# 覆盖保存
write.csv(data.frame(KeyGene = key_genes_final),
          "results/ML_key_genes_final.csv", row.names = FALSE)
save(results, imp_list, perf, ranked_imp, key_genes_final,
     file = "results/ML_results_3methods.RData")
