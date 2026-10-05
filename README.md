# GREM1-CRC multi-omics integration pipeline

Analysis scripts for the study:

**"Multi-omics prioritization of GREM1 in colorectal cancer and exploratory analysis of CAF-associated collagen communication."**

All intermediate result objects (`*.RData` / `*.rds`) were intentionally excluded from this repository. Every step is reproducible from the public datasets listed below.

## Pipeline overview

| Step | Folder / file | Content |
|---|---|---|
| 1 | `01_data_preprocessing.R` | Bulk transcriptome preprocessing and normalization (GEO, `limma`/`edgeR`, batch correction) |
| 2 | `02_MR_analysis/` | Two-sample Mendelian randomization using eQTL (`eqtl/`) and pQTL (`pqtl/`) instruments, including heterogeneity and MR-PRESSO tests |
| 3 | `03_colocalization.R` | Bayesian colocalization (`coloc`) between molecular traits and colorectal cancer risk |
| 4 | `04_single_cell_analysis/` | Single-cell analysis (Seurat/SingleR/CellChat), CAF subpopulation localization |
| 5 | `05_machine_learning/` | Machine-learning key-gene screening (caret: random forest / GLM / GBM / SVM), SHAP interpretation, consensus clustering, ssGSEA immune infiltration, survival and nomogram |
| 6 | `06_In silico single-cell perturbation/` | In-silico single-cell perturbation (scTenifoldKnk) and virtual docking screens (Python + RDKit/meeko) |

## Public data sources

- Gene expression / transcriptomes: GEO **GSE196006, GSE146009, GSE110224, GSE132465, GSE132257**
- Tissue expression reference: **GTEx v10**
- Colorectal cancer GWAS summary statistics: GWAS Catalog **GCST90013862**
- Plasma protein QTL (pQTL) instruments: Ferkingstad et al. (published cis-pQTL study)
- eQTL instruments: GTEx v10 / eQTLGen, as described in the manuscript Methods

## Requirements

### R (≥ 4.2)

Install CRAN packages and Bioconductor packages:

```r
install.packages(c(
  "TwoSampleMR", "coloc", "ieugwasr", "data.table", "dplyr", "ggplot2",
  "patchwork", "tidyr", "tidyverse", "readr", "stringr", "tibble", "purrr",
  "caret", "randomForest", "glmnet", "e1071", "fastshap", "ggbeeswarm",
  "pROC", "nnet", "survival", "survminer", "rms", "forestplot", "circlize",
  "knitr", "writexl", "reshape2", "ggrepel", "VennDiagram",
  "ConsensusClusterPlus", "GSVA", "GSEABase", "pheatmap", "arrow", "MRPRESSO"
))
BiocManager::install(c(
  "org.Hs.eg.db", "clusterProfiler", "limma", "edgeR", "GEOquery",
  "Biobase", "AnnotationDbi", "hgu133plus2.db", "sva", "GenomicRanges",
  "rtracklayer", "BiocParallel", "Matrix", "SingleR", "celldex"
))
```

GitHub-hosted packages:

```r
# Seurat (cell-level analysis)
install.packages("Seurat")              # or remotes::install_github("satijalab/seurat")
# CellChat (cell-cell communication)
remotes::install_github("jinworks/CellChat")
# scTenifoldKnk / scTenifoldNet (in-silico perturbation)
remotes::install_github("cailab-tamu/scTenifoldKnk")
remotes::install_github("cailab-tamu/scTenifoldNet")
```

### Python (≥ 3.9)

```
pandas
rdkit
meeko
```

## Run order

Follow the numeric prefixes (01 → 06). Scripts assume the public datasets listed above have been downloaded; adjust local file paths at the top of each script as needed.

## License / citation

Please cite the accompanying manuscript when using these scripts.
