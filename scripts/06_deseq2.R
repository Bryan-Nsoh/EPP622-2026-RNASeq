library(DESeq2)
library(ggplot2)
setwd("results")

samples <- read.csv("0_sample_table/samples.csv", row.names = 1)
samples$time <- factor(samples$time)

fc <- read.delim("5_counts/counts.txt", comment.char = "#", row.names = 1, check.names = FALSE)
counts <- fc[, 6:ncol(fc)]
colnames(counts) <- sub(".bam", "", basename(colnames(counts)))
rownames(counts) <- sub("gene-", "", rownames(counts))
counts <- counts[, rownames(samples)]

dds <- DESeqDataSetFromMatrix(counts, samples, design = ~ time)
dds <- DESeq(dds)
res <- results(dds, contrast = c("time", "3h", "0h"), alpha = 0.05)

png("6_deseq2/MA_plot.png")
plotMA(res)
dev.off()

res <- as.data.frame(res)
res <- res[order(res$padj), ]
sig <- subset(res, padj < 0.05 & abs(log2FoldChange) >= 1)
write.csv(res, "6_deseq2/DE_all_genes.csv")
write.csv(sig, "6_deseq2/DE_significant_padj0.05_log2FC1.csv")

vsd <- vst(dds)
saveRDS(vsd, "6_deseq2/vsd.rds")
pca <- plotPCA(vsd, intgroup = "time") + geom_text(aes(label = name), vjust = -1, show.legend = FALSE)
ggsave("6_deseq2/PCA.png", pca + xlim(-40, 40) + ylim(-10, 15), width = 8, height = 4)
