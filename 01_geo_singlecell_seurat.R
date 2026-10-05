# ============================================================
# 脚本：01_geo_singlecell_seurat.R
# 目的：下载 GSE132465/GSE132257，完成单细胞 QC、标准化、降维、聚类、UMAP 和 SingleR 注释
# 输入：GEO accession（脚本自动从 NCBI GEO 下载补充文件）
# 输出：raw/ 原始 GEO 文件；results/Seurat 对象、汇总表和图形
# 说明：脚本按 10X matrix.mtx、10X h5 和 RDS 三类常见 GEO 单细胞格式自动读取。
#       如数据集存在多个独立矩阵，会按文件夹分别建立对象并合并。
# 日期：2026-09-13
# ============================================================

required_packages <- c(
  "GEOquery", "Seurat", "SingleR", "celldex", "BiocParallel",
  "Matrix", "ggplot2", "dplyr", "readr", "writexl", "patchwork"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "缺少 R 包：", paste(missing_packages, collapse = ", "),
    "。请先运行 install_packages.R 或按 README 安装后重试。"
  )
}

suppressPackageStartupMessages({
  library(GEOquery)
  library(Seurat)
  library(SingleR)
  library(celldex)
  library(BiocParallel)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(writexl)
  library(patchwork)
})

set.seed(123)

# 路径设置：按用户指定目录组织原始文件、结果和图形。
base_dir <- "D:/Bioinfor/colorectal/single_raw"
raw_dir <- file.path(base_dir, "raw")
result_dir <- file.path(base_dir, "results")
figure_dir <- file.path(result_dir, "figures")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

accessions <- c("GSE132465", "GSE132257")

# GEO 补充文件下载到各 accession 的独立目录，便于断点重跑和溯源。
downloaded_files <- lapply(accessions, function(accession) {
  accession_dir <- file.path(raw_dir, accession)
  dir.create(accession_dir, recursive = TRUE, showWarnings = FALSE)
  GEOquery::getGEOSuppFiles(
    accession,
    makeDirectory = FALSE,
    baseDir = accession_dir,
    fetch_files = TRUE
  )
  accession_dir
})
names(downloaded_files) <- accessions

# 解包 GEO 常见 tar/zip 文件；已解包文件不会重复覆盖。
raw_archives <- list.files(raw_dir, pattern = "\\.(tar|tar\\.gz|tgz|zip)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
invisible(lapply(raw_archives, function(archive) {
  target_dir <- dirname(archive)
  if (grepl("\\.zip$", archive, ignore.case = TRUE)) {
    utils::unzip(archive, exdir = target_dir)
  } else {
    utils::untar(archive, exdir = target_dir)
  }
}))

# 将文件夹或 h5/RDS 统一转换成 Seurat 对象；这里仅保留少量必要的辅助函数。
read_matrix_object <- function(path) {
  if (dir.exists(path)) {
    matrix_file <- list.files(path, pattern = "^matrix\\.mtx(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
    if (length(matrix_file) > 0) return(Seurat::Read10X(data.dir = path))
  }
  if (grepl("\\.h5$", path, ignore.case = TRUE)) return(Seurat::Read10X_h5(path))
  if (grepl("\\.rds$", path, ignore.case = TRUE)) return(readRDS(path))
  NULL
}

make_seurat_object <- function(matrix_input, sample_id) {
  if (inherits(matrix_input, "Seurat")) {
    object <- matrix_input
    object$sample_id <- sample_id
    return(object)
  }
  if (is.list(matrix_input) && !inherits(matrix_input, "Matrix")) {
    if ("Gene Expression" %in% names(matrix_input)) matrix_input <- matrix_input[["Gene Expression"]]
    else matrix_input <- matrix_input[[1]]
  }
  Seurat::CreateSeuratObject(
    counts = matrix_input,
    project = sample_id,
    min.cells = 3,
    min.features = 200
  )
}

matrix_dirs <- list.dirs(raw_dir, recursive = TRUE, full.names = TRUE)
matrix_dirs <- matrix_dirs[vapply(matrix_dirs, function(path) {
  any(list.files(path, pattern = "^matrix\\.mtx(\\.gz)?$", full.names = FALSE, ignore.case = TRUE))
}, logical(1))]
h5_files <- list.files(raw_dir, pattern = "\\.h5$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
rds_files <- list.files(raw_dir, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
input_paths <- unique(c(matrix_dirs, h5_files, rds_files))
if (length(input_paths) == 0) {
  stop("未找到 matrix.mtx、h5 或 rds。请检查 GEO 下载是否成功及压缩包是否可解包。")
}

objects <- lapply(input_paths, function(path) {
  sample_id <- paste0(gsub("[^A-Za-z0-9]", "_", basename(dirname(path))), "_", tools::file_path_sans_ext(basename(path)))
  object <- make_seurat_object(read_matrix_object(path), sample_id)
  object$source_path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  object
})
names(objects) <- make.unique(vapply(objects, function(object) unique(object$sample_id)[1], character(1)))

# QC：线粒体比例是常用的低质量细胞指标；阈值保守，便于跨数据集复现。
objects <- lapply(objects, function(object) {
  object[["percent.mt"]] <- Seurat::PercentageFeatureSet(object, pattern = "^MT-")
  object
})
qc_before <- bind_rows(lapply(objects, function(object) {
  tibble(
    sample_id = unique(object$sample_id)[1],
    cells_before = ncol(object),
    median_features = median(object$nFeature_RNA),
    median_counts = median(object$nCount_RNA),
    median_percent_mt = median(object$percent.mt)
  )
}))
objects <- lapply(objects, function(object) {
  subset(object, subset = nFeature_RNA >= 200 & nFeature_RNA <= 7500 & percent.mt <= 20)
})

if (any(vapply(objects, ncol, integer(1)) == 0L)) stop("QC 后存在空对象，请根据 qc_before 调整阈值。")

# 合并对象后使用 SCTransform；该方法同时处理测序深度差异并用于 PCA/UMAP。
combined <- if (length(objects) == 1) objects[[1]] else Reduce(function(x, y) merge(x, y), objects)
combined <- Seurat::SCTransform(combined, vars.to.regress = "percent.mt", verbose = FALSE)
combined <- Seurat::RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- Seurat::FindNeighbors(combined, dims = 1:30, verbose = FALSE)
combined <- Seurat::FindClusters(combined, resolution = 0.5, verbose = FALSE)
combined <- Seurat::RunUMAP(combined, dims = 1:30, seed.use = 123, verbose = FALSE)

# SingleR：使用 Human Primary Cell Atlas 参考，对每个细胞给出自动标签。
reference <- celldex::HumanPrimaryCellAtlasData()
single_r <- SingleR::SingleR(
  test = Seurat::GetAssayData(combined, assay = "SCT", slot = "data"),
  ref = reference,
  labels = reference$label.main,
  BPPARAM = BiocParallel::SerialParam()
)
combined$SingleR_label <- single_r$labels
combined$SingleR_pruned_label <- single_r$pruned.labels

# 结果表：保存 QC、聚类和注释，便于后续统计或人工复核。
qc_after <- bind_rows(lapply(split(combined@meta.data, combined$sample_id), function(metadata) {
  tibble(
    sample_id = unique(metadata$sample_id),
    cells_after = nrow(metadata),
    median_features = median(metadata$nFeature_RNA),
    median_counts = median(metadata$nCount_RNA),
    median_percent_mt = median(metadata$percent.mt)
  )
}))
cluster_summary <- as.data.frame(table(combined$seurat_clusters), stringsAsFactors = FALSE) |>
  rename(cluster = Var1, cells = Freq)
annotation_summary <- as.data.frame(table(combined$SingleR_pruned_label, useNA = "ifany"), stringsAsFactors = FALSE) |>
  rename(cell_type = Var1, cells = Freq)
writexl::write_xlsx(
  list(QC_before = qc_before, QC_after = qc_after, Clusters = cluster_summary, SingleR = annotation_summary),
  file.path(result_dir, "singlecell_summary.xlsx")
)
write.csv(as.data.frame(combined@meta.data), file.path(result_dir, "cell_metadata_with_SingleR.csv"), row.names = TRUE, fileEncoding = "UTF-8")
saveRDS(combined, file.path(result_dir, "GSE132465_GSE132257_Seurat_SingleR.rds"))
writeLines(capture.output(sessionInfo()), file.path(result_dir, "sessionInfo.txt"))

# QC 图和 UMAP 图同时导出 PDF/PNG，便于论文排版和快速查看。
qc_plot <- VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.05)
ggsave(file.path(figure_dir, "QC_violin.pdf"), qc_plot, width = 12, height = 4, device = cairo_pdf)
ggsave(file.path(figure_dir, "QC_violin.png"), qc_plot, width = 12, height = 4, dpi = 300)
umap_cluster <- DimPlot(combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle("Clusters")
umap_celltype <- DimPlot(combined, reduction = "umap", group.by = "SingleR_pruned_label", label = TRUE, repel = TRUE) + ggtitle("SingleR")
umap_plot <- umap_cluster + umap_celltype
ggsave(file.path(figure_dir, "UMAP_clusters_SingleR.pdf"), umap_plot, width = 14, height = 6, device = cairo_pdf)
ggsave(file.path(figure_dir, "UMAP_clusters_SingleR.png"), umap_plot, width = 14, height = 6, dpi = 300)

stopifnot(ncol(combined) > 0, nrow(combined@meta.data) == ncol(combined), length(unique(combined$seurat_clusters)) > 1)
message("分析完成：", ncol(combined), " 个细胞，", length(unique(combined$seurat_clusters)), " 个聚类。结果位于 ", result_dir)
