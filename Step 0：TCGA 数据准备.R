setwd("D:/Bioinfor/colorectal")
rm(list = ls())
# ===== 1. 读取三个文件 =====
# 表达矩阵（HiSeqV2 无后缀，tab 分隔，第一列是基因名）
expr <- read.delim("D:/Bioinfor/colorectal/TCGA_RAW/HiSeqV2",
                   check.names = FALSE, row.names = 1)
cat("表达矩阵:", nrow(expr), "基因 ×", ncol(expr), "样本\n")
expr[1:5, 1:3]                     # 看格式
range(as.matrix(expr), na.rm = TRUE)  # 看值域（判断是否已 log2）

# 生存数据
surv <- read.delim("D:/Bioinfor/colorectal/TCGA_RAW/survival_COADREAD_survival.txt",
                   check.names = FALSE)
names(surv)                         # 应有 OS, OS.time, PFI, PFI.time 等
head(surv[, 1:5])

# 临床矩阵
clin <- read.delim("D:/Bioinfor/colorectal/TCGA_RAW/TCGA.COADREAD.sampleMap_COADREAD_clinicalMatrix",
                   check.names = FALSE)
names(clin)                         # 应有 age, gender, stage, grade 等
head(clin[, 1:5])

# ===== 2. 三表按样本 barcode 取交集 =====
# 先确认各表的样本列名
cat("表达矩阵列名示例:", colnames(expr)[1], "\n")
cat("生存第一列名:", names(surv)[1], " 示例:", surv[1, 1], "\n")
cat("临床第一列名:", names(clin)[1], " 示例:", clin[1, 1], "\n")


common_samples <- intersect(colnames(expr), 
                            intersect(surv$sample, clin$sampleID))
cat("三表交集样本数:", length(common_samples), "\n")

# 子集化对齐
#取expr中列名匹配common_samples的子集，重新赋值回expr
expr <- expr[, common_samples]
#逻辑：用match()查找common_samples在surv$sample中的位置
#效果：最终surv的样本顺序和common_samples完全对齐，只保留共有的样本
surv <- surv[match(common_samples, surv$sample), ]
clin <- clin[match(common_samples, clin$sampleID), ]

# ===== 3. 区分肿瘤 vs 正常（TCGA barcode 第 14-15 位） =====
sample_type <- substr(common_samples, 14, 15)

table(sample_type)                   # 01=肿瘤, 11=正常

tumor_samples  <- common_samples[sample_type == "01"]
normal_samples <- common_samples[sample_type == "11"]
cat("肿瘤:", length(tumor_samples), " 正常:", length(normal_samples), "\n")

# ===== 4. 保存预处理结果 =====
save(expr, surv, clin, tumor_samples, normal_samples, common_samples,
     file = "D:/Bioinfor/colorectal/TCGA/TCGA_COADREAD_processed.RData")
cat("已保存: TCGA_COADREAD_processed.RData\n")
