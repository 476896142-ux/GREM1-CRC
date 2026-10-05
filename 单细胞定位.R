# ===== 安装包（第一次运行）=====
if (!require("Seurat")) install.packages("Seurat")
if (!require("SeuratData")) remotes::install_github("satijalab/seurat-data")
if (!require("patchwork")) install.packages("patchwork")
if (!require("BiocManager")) install.packages("BiocManager")
if (!require("SingleR")) BiocManager::install("SingleR")
if (!require("celldex")) BiocManager::install("celldex")
library(Seurat)
library(SeuratData)
library(patchwork)
library(dplyr)
library(CellChat)
setwd("D:/Bioinfor/colorectal")
# ===== 下载 CRC 单细胞数据（推荐 GSE161277）=====
# 路径设置：按用户指定目录组织原始文件、结果和图形。
base_dir <- "D:/Bioinfor/colorectal/single_raw"
raw_dir <- file.path(base_dir, "raw")
result_dir <- file.path(base_dir, "results")
figure_dir <- file.path(result_dir, "figures")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

accessions <- c("GSE132465", "GSE132257")


# ===== GSE132465 =====
cat("读取 GSE132465...\n")
counts_132465 <- read.table(
  gzfile("D:/Bioinfor/colorectal/single_raw/raw/GSE132465/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz"),
  header = TRUE, row.names = 1, sep = "\t", check.names = FALSE
)
annot_132465 <- read.table(
  gzfile("D:/Bioinfor/colorectal/single_raw/raw/GSE132465/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"),
  header = TRUE, sep = "\t", check.names = FALSE, row.names = 1
)
#*构建 Seurat 对象，过滤低表达基因与低质量细胞*。
sc1 <- CreateSeuratObject(counts = counts_132465, project = "GSE132465",
                          min.cells = 3, min.features = 200)
# 加注释（确保 barcode 匹配）
common_cells <- intersect(colnames(sc1), rownames(annot_132465))
sc1 <- sc1[, common_cells]
sc1 <- AddMetaData(sc1, annot_132465[common_cells, ])
cat("GSE132465:", ncol(sc1), "细胞\n")

# ===== GSE132257（按实际文件名改）=====
cat("读取 GSE132257...\n")
counts_132257 <- read.table(
  gzfile("D:/Bioinfor/colorectal/single_raw/raw/GSE132257/GSE132257_GEO_processed_protocol_and_fresh_frozen_raw_UMI_count_matrix.txt.gz"),
  header = TRUE, row.names = 1, sep = "\t", check.names = FALSE
)

annot_132257 <- read.table(
  gzfile("D:/Bioinfor/colorectal/single_raw/raw/GSE132257/GSE132257_processed_protocol_and_fresh_frozen_cell_annotation.txt.gz"),
  header = TRUE, sep = "\t", check.names = FALSE, row.names = 1
)

sc2 <- CreateSeuratObject(counts = counts_132257, project = "GSE132257",
                          min.cells = 3, min.features = 200)
cat("GSE132257:", ncol(sc2), "细胞\n")

# ===== 合并 =====
sc <- merge(sc1, y = sc2, add.cell.ids = c("GSE132465", "GSE132257"))
cat("合并后:", ncol(sc), "细胞\n")

# ===== Seurat 对象（用 saveRDS，体积小且保留完整结构）=====
saveRDS(sc,  file.path(result_dir, "seurat_merged.rds"))      # 合并后，1.8GB
saveRDS(sc1, file.path(result_dir, "seurat_GSE132465.rds"))   # 1.5GB
saveRDS(sc2, file.path(result_dir, "seurat_GSE132257.rds"))   # 350MB

# ===== 计数矩阵（太大不用 csv，用 saveRDS）=====
saveRDS(counts_132465, file.path(result_dir, "counts_GSE132465.rds"))
saveRDS(counts_132257, file.path(result_dir, "counts_GSE132257.rds"))

# ===== 细胞注释（小，用 csv 方便查看）=====
write.csv(annot_132465, file.path(result_dir, "annotation_GSE132465.csv"))
write.csv(annot_132257, file.path(result_dir, "annotation_GSE132257.csv"))

cat("全部保存完成！\n")
cat("目录:", result_dir, "\n")
list.files(result_dir)
# ================================================================
# 2. QC + 标准化 + 降维聚类
# ================================================================
#抽样
set.seed(123)
sc <- subset(sc, cells = sample(colnames(sc), 25000))  # 抽 2.5 万细胞
cat("抽样后:", ncol(sc), "细胞 ×", nrow(sc), "基因\n")
#统计每个细胞中 ** 线粒体基因（MT - 开头）** 的 reads 占总 reads 的百分比
sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^MT-")
#基于 Seurat 对象sc做细胞过滤，筛选合格细胞;筛选条件：基因数nFeature_RNA：200 ＜ 检测基因数 ＜ 6000
#线粒体基因占比percent.mt ＜ 20%
#目的：去除低质量、破损、双细胞，保留可靠单细胞用于下游分析
sc <- subset(sc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 20)
cat("QC后:", ncol(sc), "细胞\n")
#LogNormalize（对数归一化）,对单细胞对象sc执行标准化
sc <- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)
#Seurat 筛选单细胞可变基因
sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)
#对指定基因表达量标准化 + 中心化缩放（均值归零、方差归一）,用于后续 PCA,只对高变基因做 ScaleData
sc <- ScaleData(sc, features = VariableFeatures(sc))
#对 Seurat 对象sc执行 PCA 主成分分析
sc <- RunPCA(sc, features = VariableFeatures(sc), verbose = FALSE)
#基于 PCA 降维结果，计算细胞间 K 近邻，构建细胞邻接图，用于后续聚类
sc <- FindNeighbors(sc, dims = 1:30)
#Seurat 包细胞聚类
sc <- FindClusters(sc, resolution = 0.5)
#对单细胞数据执行 UMAP 降维，结果存入 sc 对象，用于可视化细胞分群
sc <- RunUMAP(sc, dims = 1:30, verbose = FALSE)

# 保存聚类 UMAP
p1 <- DimPlot(sc, reduction = "umap", label = TRUE) + ggtitle("Seurat Clusters")
ggsave(file.path(figure_dir, "01_umap_clusters.png"), p1, width = 8, height = 6)

# ================================================================
# 3. 细胞类型注释
# ================================================================
# 优先用 GSE132465 自带的 celltype 列（如果有）
if ("cell_type" %in% colnames(sc@meta.data)) {
  sc$celltype <- sc$cell_type
} else if ("Cell_type" %in% colnames(sc@meta.data)) {
  sc$celltype <- sc$Cell_type
} else {
  # 没有自带注释就用 SingleR 自动注释
  library(SingleR)
  library(celldex)
  ref <- celldex::HumanPrimaryCellAtlasData()
  sc_data <- GetAssayData(sc, slot = "data")
  pred <- SingleR(test = sc_data, ref = ref, labels = ref$label.main,
                  clusters = sc$seurat_clusters)
  sc$celltype <- pred$labels[match(sc$seurat_clusters, rownames(pred))]
}

p2 <- DimPlot(sc, reduction = "umap", group.by = "celltype",
              label = TRUE, repel = TRUE) + ggtitle("Cell Types")
ggsave(file.path(figure_dir, "02_umap_celltype.png"), p2, width = 10, height = 7)

# 保存细胞类型统计
celltype_stats <- as.data.frame(table(sc$celltype))
colnames(celltype_stats) <- c("CellType", "n_cells")
write.csv(celltype_stats, file.path(result_dir, "celltype_stats.csv"), row.names = FALSE)
print(celltype_stats)

# ================================================================
