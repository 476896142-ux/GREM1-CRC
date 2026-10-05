Sys.setenv(OPENGWAS_JWT = "eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiI0NzY4OTYxNDJAcXEuY29tIiwiaWF0IjoxNzg4MDgzMjg2LCJleHAiOjE3ODkyOTI4ODZ9.Fyt4awhW1fazfRlOqfqNChCzdqfnclBHHt5z8qCQ32ZSbFcoYC2R1OExZSFmbVqHVMAjgphVsmlU6HjVKUzQ9MTHoNZgPs8Bko0qjJr1MCSaudbIxBU89kv0wKsq9VstNdS7efXEGais48v4nToTHks3BUwWPdNSKKD1AbY9tcDjZ8Tp73kyA6TTb8OfrGtOYuBDfImQJhz4HpzEWiGi-fnUuBT0Z7PI9lQ2v4YFr9jaf4HHUwyynDZlRh9qlpjrplzGTFK65cngWzvCbeA87J74YT08nTX1BMGpl7PidXXzRzySAIEzPPbLf3XWUhbY2wv4uLJVEY0VOH_llefSew")

suppressPackageStartupMessages({
  library(arrow); library(data.table); library(dplyr); library(tidyr); library(stringr)
  library(MungeSumstats); library(TwoSampleMR);library(MendelR);library(QTLMR)
})
set.seed(123)
project_dir <- Sys.getenv("CRC_MR_DIR", "D:/Bioinfor/colorectal/MR")
output_dir <- file.path(project_dir, "Mr_colorectal_result")
sig_file <- file.path(project_dir, "Colon_Sigmoid.v10.eQTLs.signif_pairs.parquet")
tra_file <- file.path(project_dir, "Colon_Transverse.v10.eQTLs.signif_pairs.parquet")
ld_ref <- file.path(project_dir, "g1000_eur", "g1000_eur")
plink_path <- file.path(project_dir, "LD", "plink_win64", "plink.exe")
p_thresh <- 5e-6; 
r2_thresh <- 0.1; 
kb_window <- 10000; 
f_thresh <- 10
message("[1/6] 检查输入路径和 PLINK reference")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
required_files <- c(sig_file, tra_file, plink_path, paste0(ld_ref, c(".bed", ".bim", ".fam")))
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) stop("缺少输入文件：\n", paste(missing_files, collapse = "\n"))

message("[2/6] 读取并合并两个组织")
required_cols <- c("gene_id", "variant_id", "pval_nominal", "af", "slope", "slope_se","ma_count")

eqtl_data_sig <- arrow::read_parquet(sig_file, col_select = required_cols) |> mutate(tissue = "Colon_Sigmoid") |> as.data.table()

eqtl_data_tra <- arrow::read_parquet(tra_file, col_select = required_cols) |> mutate(tissue = "Colon_Transverse") |> as.data.table()

missing_cols <- setdiff(required_cols, union(names(eqtl_data_sig), names(eqtl_data_tra)))
if (length(missing_cols)) stop("eQTL 缺少必要列：", paste(missing_cols, collapse = ", "))
eqtl_all <- rbindlist(list(eqtl_data_sig, eqtl_data_tra), use.names = TRUE, fill = TRUE)
eqtl_all[, n_tissue := uniqueN(tissue), by = .(gene_id, variant_id)]
eqtl_dedup <- eqtl_all[order(pval_nominal), .SD[1L], by = .(gene_id, variant_id)]

message("[3/6] 整理字段并转换 rsID")
View(harmonised_data)
eqtl_data <- eqtl_dedup %>%
  separate(variant_id, into = c("CHR", "BP", "A1", "A2", "Build"), sep = "_",
           remove = FALSE, extra = "merge", fill = "right") %>%
  mutate(CHR = str_remove(as.character(CHR), "^chr"), BP = as.integer(BP),
         A1 = toupper(A1), A2 = toupper(A2), EAF = as.numeric(af),
         BETA = as.numeric(slope), SE = as.numeric(slope_se), P = as.numeric(pval_nominal),
         ma_count = as.numeric(ma_count),
         maf = pmin(EAF, 1 - EAF),
         samplesize = ma_count / (2 * maf))%>%
  filter(!is.na(CHR), !is.na(BP), !is.na(A1), !is.na(A2), is.finite(BETA),
         is.finite(SE), SE > 0, is.finite(P), P > 0, P <= 1,
         is.finite(ma_count), ma_count >= 0, maf > 0, is.finite(samplesize)
  )
View(eqtl_data)
make_key <- function(chr, bp, a1, a2) paste0(tolower(sub("^chr", "", chr)), ":", bp, ":",
                                             toupper(a1), ":", toupper(a2))
eqtl_data <- mutate(eqtl_data, join_key = make_key(CHR, BP, A1, A2))

formatted <- MungeSumstats::format_sumstats(
  as.data.frame(select(eqtl_data, all_of(c("CHR", "BP", "A1", "A2", "BETA", "P","samplesize")))),
  ref_genome = "GRCh38",
  nThread = max(1L, parallel::detectCores() - 1L), dbSNP = 155, return_data = TRUE,
  on_ref_genome = TRUE) |>
  mutate(join_key = make_key(CHR, BP, A1, A2)) |>
  left_join(select(eqtl_data, join_key, gene_id, tissue, n_tissue, EAF, SE), by = "join_key") |>
  filter(!is.na(gene_id))

if (!nrow(formatted)) stop("rsID 转换后无法回接 gene_id，请检查 variant_id 格式")

exposure_dat <- TwoSampleMR::format_data(
  as.data.frame(formatted),
  type = "exposure", snp_col = "SNP", beta_col = "BETA", se_col = "SE",
  effect_allele_col = "A2", other_allele_col = "A1", pval_col = "P", eaf_col = "EAF",samplesize_col = "N",
  phenotype_col = "gene_id") |>
  mutate(id.exposure = paste0("Colorectal eQTL: ", exposure), F_stat = (beta.exposure / se.exposure)^2) |>
  filter(pval.exposure <= p_thresh, is.finite(F_stat), F_stat >= f_thresh, !is.na(SNP), SNP != "") |>
  distinct(id.exposure, SNP, .keep_all = TRUE)

if (!nrow(exposure_dat)) stop("显著性和 F 统计量筛选后没有工具变量")

write.table(exposure_dat, file.path(output_dir, "Colorectal_exposure_preclump.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

clump_res <- ieugwasr::ld_clump(
  dat = exposure_dat |>
    rename(rsid = SNP, pval = pval.exposure),    # ← 改这两列名
  clump_kb  = kb_window,
  clump_r2  = r2_thresh,
  clump_p   = p_thresh,
  bfile     = ld_ref,
  plink_bin = plink_path
)
View(clump_res)
# clump 完改回 TwoSampleMR 标准列名
exp_dat_clump <- clump_res |>
  rename(SNP = rsid, pval.exposure = pval)
message("Clumping 完成: ", n_distinct(exp_dat_clump$exposure), " 个基因, ",
        nrow(exp_dat_clump), " 个 SNP")
# ========== 3. 弱工具筛选（F > 10） ==========
exp_dat_final <- exp_dat_clump %>%
  filter(F_stat > f_thresh)
# ========== 4. 结果汇总与导出 ==========
exp_dat_final <- exp_dat_final %>%
  left_join(
    formatted %>% select(SNP, CHR) %>% distinct(SNP, .keep_all = TRUE),
    by = "SNP"
  ) %>%
  rename(chr.exposure = CHR)
readr::write_csv(exp_dat_final, "IVs_selected_colorectal.csv")

View(exp_dat_final)
message("Clumping 完成: ", n_distinct(exp_dat_clump$exposure), " 个基因, ",
        nrow(exp_dat_clump), " 个 SNP")
exp_dat<- readr::read_csv("IVs_selected_colorectal.csv", show_col_types = FALSE)
outcome_raw <- extract_outcome_data(
  snps     = exp_dat$SNP,        # 字符向量，API 路径按 rsID 查，不受 build 影响
  outcomes = "ebi-a-GCST90013862",
  proxies  = F
) |>
  mutate(id.outcome = "ebi-a-GCST90013862",
         outcome = "Colorectal cancer",
         samplesize.outcome = 407746,
         ncase.outcome = 4161,
         ncontrol.outcome = 403585)
View(outcome_raw)

harmonised_data <- harmonise_data(exp_dat, outcome_raw) |>
  filter(mr_keep)

# 增加file参数
save(eqtl_all, eqtl_dedup, eqtl_data, exp_dat, exp_dat_clump, exp_dat_final, outcome_raw, harmonised_data, file = "eqtl_result.RData")
