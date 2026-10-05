# ================================================================
# C阶段：共定位分析
# pQTL 基因用 p1e-05 完整数据，eQTL 基因用 GTEx v10 结肠 eQTL
# ================================================================
setwd("D:/Bioinfor/colorectal")
library(coloc)
library(dplyr)
library(TwoSampleMR)
library(org.Hs.eg.db)
library(ieugwasr)
# ===== 基因分组 =====
pqtl_genes <- c("GREM1","HHIP","LAMC2","NTRK3","ANGPTL1","KLB")
eqtl_genes <- c("CSE1L","EPB41L3","EPHB1","WDR72")
WINDOW <- 500000

# ================================================================
# 第一部分：pQTL 基因共定位
# ================================================================


pqtl_05 <- readr::read_csv("D:/Bioinfor/colorectal/MR/pQTL/exposure_data.csv",
                           show_col_types = FALSE)
pqtl_coloc_results <- list()

for (gene in pqtl_genes) {
  cat("\n--- pQTL:", gene, "---\n")
  
  gene_pqtl <- pqtl_05 %>%
    filter(exposure == gene) %>%
    filter(!is.na(SNP), SNP != "", !duplicated(SNP),
           !is.na(beta.exposure), !is.na(se.exposure), se.exposure > 0)
  
  cat("SNP 数:", nrow(gene_pqtl), "\n")
  if (nrow(gene_pqtl) < 10) { cat("SNP 太少，跳过\n"); next }
  
  gwas_region <- tryCatch({
    extract_outcome_data(snps = gene_pqtl$SNP,
                         outcomes = "ebi-a-GCST90013862", proxies = FALSE)
  }, error = function(e) { cat("GWAS 提取失败:", e$message, "\n"); NULL })
  
  if (is.null(gwas_region) || nrow(gwas_region) == 0) next
  
  # 关键改动：as.data.frame() 转换TwoSampleMR对象
  merged <- gene_pqtl %>%
    inner_join(
      as.data.frame(gwas_region) %>% dplyr::select(SNP, beta.outcome, se.outcome),
      by = "SNP"
    ) %>%
    filter(!is.na(beta.outcome), !is.na(se.outcome), se.outcome > 0) %>%
    distinct(SNP, .keep_all = TRUE)
  
  cat("对齐后 SNP 数:", nrow(merged), "\n")
  if (nrow(merged) < 8) next
  
  pqtl_N <- round(mean(merged$samplesize.exposure, na.rm = TRUE))
  
  res <- coloc.abf(
    dataset1 = list(snp = merged$SNP, beta = merged$beta.exposure,
                    varbeta = merged$se.exposure^2, type = "quant",
                    N = pqtl_N,MAF = merged$eaf.exposure),
    dataset2 = list(snp = merged$SNP, beta = merged$beta.outcome,
                    varbeta = merged$se.outcome^2, type = "cc",
                    N = 407746, s = 4161/407746)
  )
  
  pqtl_coloc_results[[gene]] <- res$summary
  cat("PP.H4 =", round(res$summary["PP.H4.abf"], 3),
      "  PP.H3 =", round(res$summary["PP.H3.abf"], 3), "\n")
}

save(pqtl_coloc_results, file = "pqtl_coloc_result.RData")
load("pqtl_coloc_result.RData")


# 第二部分：eQTL 基因共定位
# ================================================================
cat("\n========== eQTL 共定位 ==========\n")
if (!require("BiocManager")) install.packages("BiocManager")
if (!require("rtracklayer")) BiocManager::install("rtracklayer")
if (!require("GenomicRanges")) BiocManager::install("GenomicRanges")
library(rtracklayer)
library(GenomicRanges)
library(dplyr)
library(ieugwasr)
library(coloc)
chain_file <- "D:/Bioinfor/colorectal/MR/hg38ToHg19.over.chain"
chain <- import.chain(chain_file)
cat("chain 文件加载完成\n")

load("D:/Bioinfor/colorectal/MR/eQTL/eqtl_result.RData")

eqtl_gene_map <- AnnotationDbi::select(org.Hs.eg.db,
                                       keys = eqtl_genes, keytype = "SYMBOL",
                                       columns = c("ENSEMBL","SYMBOL")) %>%
  filter(!is.na(ENSEMBL)) %>%
  mutate(ENSEMBL = sub("\\..*$", "", ENSEMBL))

eqtl_coloc_results <- list()

for (i in 1:nrow(eqtl_gene_map)) {
  gene_symbol <- eqtl_gene_map$SYMBOL[i]
  gene_ensembl <- eqtl_gene_map$ENSEMBL[i]
  cat("\n--- eQTL:", gene_symbol, "(", gene_ensembl, ") ---\n")
  
  gene_eqtl <- eqtl_data %>%
    mutate(gene_clean = sub("\\..*$", "", gene_id)) %>%
    filter(gene_clean == gene_ensembl,
           !is.na(BETA), !is.na(SE), SE > 0, !is.na(EAF)) %>%
    filter(!duplicated(paste0(CHR, ":", BP)))
  
  cat("GRCh38 SNP 数:", nrow(gene_eqtl), "\n")
  if (nrow(gene_eqtl) < 8) { cat("SNP 太少，跳过\n"); next }
  
  # 构造 GRanges（GRCh38）并 liftOver 到 GRCh37
  gr_38 <- GRanges(
    seqnames = paste0("chr", gene_eqtl$CHR),
    ranges = IRanges(start = gene_eqtl$BP, end = gene_eqtl$BP),
    index = 1:nrow(gene_eqtl)
  )
  gr_37_list <- liftOver(gr_38, chain)
  keep <- elementNROWS(gr_37_list) == 1
  gr_37_mapped <- unlist(gr_37_list[keep])
  
  gene_eqtl_mapped <- gene_eqtl[keep, ]
  gene_eqtl_mapped$CHR37 <- sub("^chr", "", as.character(seqnames(gr_37_mapped)))
  gene_eqtl_mapped$BP37 <- start(gr_37_mapped)
  gene_eqtl_mapped$chrpos37 <- paste0(gene_eqtl_mapped$CHR37, ":", gene_eqtl_mapped$BP37)
  
  cat("liftOver 成功:", nrow(gene_eqtl_mapped), "个 SNP\n")
  if (nrow(gene_eqtl_mapped) < 8) next
  
  # 提取 GWAS（分批）
  chrpos_list <- gene_eqtl_mapped$chrpos37
  gwas_list <- list()
  for (chunk in split(chrpos_list, ceiling(seq_along(chrpos_list)/100))) {
    gwas_chunk <- tryCatch({
      ieugwasr::associations(variants = chunk, id = "ebi-a-GCST90013862")
    }, error = function(e) { cat("GWAS 提取失败:", e$message, "\n"); NULL })
    if (!is.null(gwas_chunk) && nrow(gwas_chunk) > 0) {
      gwas_list[[length(gwas_list)+1]] <- gwas_chunk
    }
  }
  gwas_all <- do.call(rbind, gwas_list)
  if (is.null(gwas_all) || nrow(gwas_all) == 0) next
  
  names(gwas_all) <- tolower(names(gwas_all))
  pos_col <- if ("position" %in% names(gwas_all)) "position" else "pos"
  gwas_all <- gwas_all %>%
    mutate(chrpos37 = paste0(chr, ":", .data[[pos_col]])) %>%
    filter(!duplicated(chrpos37))
  
  merged <- gene_eqtl_mapped %>%
    inner_join(gwas_all %>% dplyr::select(chrpos37, beta, se), by = "chrpos37") %>%
    filter(!is.na(beta), !is.na(se), se > 0) %>%
    distinct(chrpos37, .keep_all = TRUE)
  
  cat("对齐后 SNP 数:", nrow(merged), "\n")
  if (nrow(merged) < 8) next
  
  eqtl_N <- round(mean(merged$samplesize, na.rm = TRUE))
  
  res <- coloc.abf(
    dataset1 = list(snp = merged$chrpos37, beta = merged$BETA,
                    varbeta = merged$SE^2, type = "quant",
                    N = eqtl_N, MAF = merged$EAF),
    dataset2 = list(snp = merged$chrpos37, beta = merged$beta,
                    varbeta = merged$se^2, type = "cc",
                    N = 407746, s = 4161/407746)
  )
  
  eqtl_coloc_results[[gene_symbol]] <- res$summary
  cat("PP.H4 =", round(res$summary["PP.H4.abf"], 3),
      "  PP.H3 =", round(res$summary["PP.H3.abf"], 3), "\n")
}

save(eqtl_coloc_results, file = "eqtl_coloc_result.RData")
load("eqtl_coloc_result.RData")
# ================================================================
# 第三部分：合并结果
# ================================================================
all_coloc <- c(pqtl_coloc_results, eqtl_coloc_results)

if (length(all_coloc) > 0) {
  coloc_df <- do.call(rbind, all_coloc) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Gene") %>%
    mutate(Source = ifelse(Gene %in% pqtl_genes, "pQTL", "eQTL")) %>%
    arrange(desc(PP.H4.abf))
  
  cat("\n===== 共定位结果汇总 =====\n")
  print(as.data.frame(coloc_df[, c("Gene","Source","PP.H3.abf","PP.H4.abf")]),
        row.names = FALSE)
  
  write.csv(coloc_df, "D:/Bioinfor/colorectal/MR/results/coloc_results.csv",
            row.names = FALSE)
  
  cat("\n强共定位（PP.H4>0.8）:",
      paste(coloc_df$Gene[coloc_df$PP.H4.abf > 0.8], collapse = ", "), "\n")
  cat("结果已保存: MR/results/coloc_results.csv\n")
} else {
  cat("没有基因完成共定位\n")
}

