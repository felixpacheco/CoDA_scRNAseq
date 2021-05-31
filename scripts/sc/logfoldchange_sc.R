# libraries
library("tidyverse")
library("xtable")
# load the 3 log-fold change
deseq <- read_csv("/Users/laurasansc/Desktop/BioinformaticsSystemsBiology/CoDA_special_course/SC_deseq2_DE.csv.gz")
edge <- read_csv("/Users/laurasansc/Desktop/BioinformaticsSystemsBiology/CoDA_special_course/SC_edgeR_DE.csv.gz")
aldex <- read.csv("results/bulk/aldex2_DE.csv")

deseq <- deseq[, c(1, 2, 3, 4, 5, 6, 7)]
edge <- edge[, c(1, 2, 3, 4, 5, 6, 7)]
# -------------------- START  VOLCANO PLOTS ---------------------------------
png(file = "volcanoplots.png", height = 1500, width = 4000, res = 300)
par(mfrow = c(1, 2))

# Volcano plot desqe2
with(deseq, plot(log2FoldChange, -log10(pvalue), pch = 1, main = "Volcano plot: DESeq2", xlab = "Log2FoldChange", ylab = "Log10(pvalue)", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5))
# with(subset(deseq, padj>0.5), points(log2FoldChange, -log10(pvalue), pch=20, col="gray"))
with(subset(deseq, log2FoldChange < -1 & padj < 0.5), points(log2FoldChange, -log10(pvalue), pch = 1, col = "red"))
with(subset(deseq, log2FoldChange > 1 & padj < 0.5), points(log2FoldChange, -log10(pvalue), pch = 1, col = "blue"))
abline(h = 1, col = "blue", lty = 2, lwd = 1)
abline(v = c(-1, 1), col = "blue", lty = 2, lwd = 1)

# Volcano plot EdgeR
with(edge, plot(logFC, -log10(PValue), pch = 1, main = "Volcano plot: EdgeR", xlab = "Log2FoldChange", ylab = "Log10(pvalue)", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5))
# with(subset(edge, FDR>0.5), points(logFC, -log10(PValue), pch=20, col="gray"))
with(subset(edge, logFC < -1 & FDR < 0.5), points(logFC, -log10(PValue), pch = 1, col = "red"))
with(subset(edge, logFC > 1 & FDR < 0.5), points(logFC, -log10(PValue), pch = 1, col = "blue"))
abline(h = 1, col = "blue", lty = 2, lwd = 1)
abline(v = c(-1, 1), col = "blue", lty = 2, lwd = 1)


# Volcano plot ALDEx2
with(aldex, plot(diff.btw, -log10(we.ep), pch = 1, main = "Volcano plot: ALDEx2", xlab = "Log2FoldChange", ylab = "Log10(pvalue)", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5))
# with(subset(aldex, we.eBH>0.5), points(diff.btw, -log10(we.ep), pch=20, col="gray"))
with(subset(aldex, diff.btw < -1 & we.eBH < 0.5), points(diff.btw, -log10(we.ep), pch = 1, col = "red"))
with(subset(aldex, diff.btw > 1 & we.eBH < 0.5), points(diff.btw, -log10(we.ep), pch = 1, col = "blue"))
abline(h = 1, col = "blue", lty = 2, lwd = 1)
abline(v = c(-1, 1), col = "blue", lty = 2, lwd = 1)

# dev.off()
dev.off()

# tables of subsets
deseq_up <- subset(deseq, log2FoldChange > 1 & padj < 0.5)
edge_up <- subset(edge, logFC > 1 & FDR < 0.5)
aldex_up <- subset(aldex, diff.btw > 1 & we.eBH < 0.5)

deseq_up_c <- data.frame("X" = deseq_up$X, "DESeq2" = deseq_up$log2FoldChange)
edge_up_c <- data.frame("X" = edge_up$X, "EdgeR" = edge_up$logFC)
aldex_up_c <- data.frame("X" = aldex_up$X, "ALDEx2" = aldex_up$diff.win)

d <- merge(deseq_up_c, edge_up_c, by = c("X"), all = TRUE)
d <- merge(d, aldex_up_c, by = c("X"), all = TRUE)

boxplot(d[, -1])
# to get the number of upregulated genes
library(reshape2)
library(ggthemes)

d <- melt(d)

png(file = "boxplot_logfoldchange.png", height = 1200, width =2500, res = 300)
ggplot(data = d, aes(x = variable, y = value)) +
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.1) +
  theme_base() +
  labs(x = "DE Package", y = "Log2FoldChange") +
  annotate("text", x = 1, y = 11, label = "1490") +
  annotate("text", x = 2, y = 11, label = "1061") +
  annotate("text", x = 3, y = 11, label = "17")
dev.off()  

