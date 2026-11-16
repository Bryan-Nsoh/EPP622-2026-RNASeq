library(DESeq2)
library(pheatmap)
setwd("results")

samples <- read.csv("0_sample_table/samples.csv", row.names = 1)
vsd <- readRDS("6_deseq2/vsd.rds")
sig <- read.csv("6_deseq2/DE_significant_padj0.05_log2FC1.csv", row.names = 1)
top50 <- head(rownames(sig), 50)

pheatmap(assay(vsd)[top50, ], scale = "row", annotation_col = samples["time"], fontsize_row = 6,
         filename = "7_heatmap/top50_DE_genes_heatmap.png", width = 7, height = 10)
