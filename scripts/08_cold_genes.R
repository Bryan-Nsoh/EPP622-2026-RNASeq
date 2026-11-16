setwd("results")

res <- read.csv("6_deseq2/DE_all_genes.csv", row.names = 1)

cold <- c(CBF1 = "AT4G25490", CBF2 = "AT4G25470", CBF3 = "AT4G25480",
          COR15A = "AT2G42540", COR15B = "AT2G42530", COR47 = "AT1G20440",
          COR78 = "AT5G52310", COR6.6 = "AT5G15970", COR413 = "AT2G15970")

cold_genes <- res[cold, c("baseMean", "log2FoldChange", "padj")]
cold_genes$name <- names(cold)
write.csv(cold_genes, "8_cold_genes/COR_CBF_log2FC.csv")
