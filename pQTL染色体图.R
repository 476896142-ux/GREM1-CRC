library(circlize)

circos.clear()

cyto_path <- system.file(package = "circlize", "extdata", "cytoBand.txt")
cyto_band <- read.cytoband(cyto_path)$df

mr_dir <- "D:/Bioinfor/colorectal/MR"
result_dir <- file.path(mr_dir, "results", "GCST90013862")

if (!dir.exists(result_dir)) {
  dir.create(result_dir, recursive = TRUE)
}

snp_df <- data.frame(
  chr    = c("chr15","chr4","chr1","chr10","chr15","chr1","chr4","chr4"),
  start  = c(32994756,145659064,182968079,7648190,88794859,178497186,39504269,39413780),
  end    = c(32994756,145659064,182968079,7648190,88794859,178497186,39504269,39413780),
  value  = rep(1, 8)
)

out_file <- file.path(result_dir, "circos_gene_snp_circle.png")

png(out_file, width = 2400, height = 2400, res = 300)
par(mar = c(0, 0, 0, 0))

circos.initializeWithIdeogram(
  cyto_band,
  plotType = c("ideogram", "labels")
)

circos.genomicTrack(
  snp_df[, c("chr", "start", "end", "value")],
  ylim = c(0, 1),
  track.height = 0.08,
  cell.padding = c(0, 0, 0, 0),
  panel.fun = function(region, value, ...) {
    circos.genomicPoints(
      region, value,
      pch = 16,
      cex = 1.4,
      col = "#E63946"
    )
  }
)

dev.off()

cat("图片已保存到：", out_file, "\n")

# 看完图之后，手动单独运行 circos.clear() 释放画布
#################线性染色体图#####################
library(ggplot2)

snp_df <- data.frame(
  chr    = c("chr15","chr4","chr1","chr10","chr15","chr1","chr4","chr4"),
  snp_pos= c(32994756,145659064,182968079,7648190,88794859,178497186,39504269,39413780),
  gene_pos = c(92855861, 1466812, 134601606,126889400, 88712000,179718500, 34442000, 37927500),
  Gene   = c("GREM1","HHIP","LAMC2","ITIH5","NTRK3","ANGPTL1","UGDH","KLB")
)

# 去掉chr前缀，转为数字
snp_df$chr_num <- as.integer(gsub("chr","",snp_df$chr))

ggplot(snp_df) +
  geom_vline(aes(xintercept = snp_pos), linetype = 2, color="#E63946", alpha=0.7) +
  geom_point(aes(x = gene_pos, y = chr_num), size=4, color="#2c3e50") +
  geom_text(aes(x = gene_pos, y = chr_num, label=Gene), nudge_y = 0.35, size=3.5) +
  scale_y_reverse(breaks = 1:22) +
  labs(x = "Genomic position (hg19, bp)", y = "Chromosome") +
  theme_bw() +
  theme(panel.grid = element_blank())

ggsave(
  file.path(result_dir, "linear_chromosome_gene_snp.png"),
  width = 14, 
  height = 10, 
  dpi = 300
)
