# 拼接 Fig3_4_3 (A) + Fig3_4_4 (B)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(ggbeeswarm)

out_dir <- "results"

# ===== 重新读图3（融合排名） =====
# 重新加载数据
rm(list = ls())
setwd("D:/Bioinfor/colorectal")
load("TCGA/TCGA_COADREAD_processed.RData")

genes_info <- read.csv("MR/results/all_genes_with_DEG_info.csv")
candidate_genes <- genes_info$GeneSymbol
candidate_genes[candidate_genes == "MS4A8"] <- "MS4A8B"
candidate_genes <- candidate_genes[candidate_genes != "GPER1"]
candidate_genes <- intersect(candidate_genes, rownames(expr))
X <- as.matrix(t(expr[candidate_genes, , drop = FALSE]))

# 重新跑ML重要性和融合排名（简化版，用之前的结果）
# 直接读之前保存的排名结果
ranked_imp <- read.csv("results/ML_key_genes_ranked.csv")
ppi <- read.csv("MR/results/hub_genes_39nodes.csv")
ppi_degree <- ppi %>% select(GeneSymbol, Degree) %>% rename(Gene = GeneSymbol)

w_ml  <- 1
w_ppi <- 2.5
combined <- ranked_imp %>%
  left_join(ppi_degree, by = "Gene") %>%
  mutate(Degree = ifelse(is.na(Degree), 0, Degree),
         ML_z   = as.numeric(scale(MeanImportance)),
         PPI_z  = as.numeric(scale(Degree)),
         Combined = w_ml * ML_z + w_ppi * PPI_z) %>%
  arrange(desc(Combined))

key10 <- combined %>% head(10) %>% pull(Gene)

# ===== Panel A: 融合排名 =====
combined_plot <- combined %>% head(15) %>%
  mutate(Gene = factor(Gene, levels = rev(Gene)),
         is_key = ifelse(Gene %in% key10, "Key gene (Top10)", "Other"))
pA <- ggplot(combined_plot, aes(x = Combined, y = Gene, fill = is_key)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Key gene (Top10)" = "#E63946",
                               "Other" = "#A8DADC"), name = NULL) +
  geom_text(aes(label = sprintf("%.2f", Combined)),
            hjust = -0.15, size = 3.5) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
  labs(title = "A", x = "Combined z-score (ML + PPI)", y = NULL) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 18),
        panel.grid.major.y = element_blank(),
        legend.position = c(0.7, 0.08),
        legend.background = element_rect(fill = "white", color = NA))

# ===== Panel B: SHAP蜂群图（重新跑SHAP） =====
library(caret)
library(randomForest)
library(fastshap)

# 重新跑RF
ctrl <- trainControl(method = "cv", number = 5, classProbs = TRUE,
                     summaryFunction = twoClassSummary, savePredictions = "final")
sample_type <- substr(rownames(X), 14, 15)
y <- factor(ifelse(sample_type == "01", "Tumor", "Normal"),
            levels = c("Normal", "Tumor"))
set.seed(123)
rf_fit <- train(X, y, method = "rf", trControl = ctrl, metric = "ROC",
                tuneGrid = data.frame(mtry = 4))
rf_final <- rf_fit$finalModel

set.seed(456)
idx <- sample(seq_len(nrow(X)), min(300, nrow(X)))
X_explain <- X[idx, , drop = FALSE]
pred_wrapper <- function(object, newdata) {
  predict(object, newdata = newdata, type = "prob")[, "Tumor"]
}
set.seed(789)
shap_res <- fastshap::explain(object = rf_final, X = X_explain,
                              pred_wrapper = pred_wrapper, nsim = 10)
shap_long <- shap_res %>% as.data.frame() %>%
  mutate(sample = row_number()) %>%
  pivot_longer(-sample, names_to = "Gene", values_to = "SHAP")
expr_long <- as.data.frame(X_explain) %>%
  mutate(sample = row_number()) %>%
  pivot_longer(-sample, names_to = "Gene", values_to = "Expression")
shap_plot_df <- shap_long %>% left_join(expr_long, by = c("sample", "Gene")) %>%
  group_by(Gene) %>% mutate(mean_abs_shap = mean(abs(SHAP))) %>% ungroup() %>%
  mutate(Gene = factor(Gene, levels = unique(Gene[order(mean_abs_shap)])))

pB <- ggplot(shap_plot_df, aes(x = SHAP, y = Gene, colour = Expression)) +
  geom_quasirandom(size = 1.6, width = 0.4, alpha = 0.8) +
  scale_color_gradientn(colours = c("#2166AC", "white", "#B2182B"),
                        name = "Expression") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  labs(title = "B", x = "SHAP value (impact on model output)", y = NULL) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 18),
        panel.grid.major.y = element_blank())

# ===== 横向拼接 =====
p_combined <- pA + pB + plot_layout(ncol = 2, widths = c(1, 1.2))

ggsave(file.path(out_dir, "Fig3_4_AB_combined_SHAP.png"),
       p_combined, width = 16, height = 9, dpi = 300)
cat("拼接完成: results/Fig3_4_AB_combined_SHAP.png\n")
