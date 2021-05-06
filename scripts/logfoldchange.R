# libraries
library("tidyverse")
library("xtable")
# load the 3 log-fold change
deseq <- read.csv(unz("results/bulk/deseq2_DE.csv.zip", "deseq2_DE.csv"))
edge <- read.csv(unz("results/bulk/edgeR_DE.csv.zip", "edgeR_DE.csv"))
aldex <- read.csv("results/bulk/aldex2_DE.csv")

deseq <- deseq[, c(1, 2, 3, 4, 5, 6, 7)]
edge <-edge[, c(1, 2, 3, 4, 5, 6, 7)]
# -------------------- START  VOLCANO PLOTS ---------------------------------
png(file = "volcanoplots.png", height = 1500, width = 4500, res = 300)
par(mfrow=c(1,3))

# Volcano plot desqe2
with(deseq, plot(log2FoldChange, -log10(pvalue),pch=1, main="Volcano plot: DESeq2",xlab="Log2FoldChange", ylab="Log10(pvalue)"))
#with(subset(deseq, padj>0.5), points(log2FoldChange, -log10(pvalue), pch=20, col="gray"))
with(subset(deseq, log2FoldChange < -1 & padj<0.5), points(log2FoldChange, -log10(pvalue), pch=1, col="red"))
with(subset(deseq, log2FoldChange > 1 & padj<0.5), points(log2FoldChange, -log10(pvalue), pch=1, col="green"))
abline(h = 1, col = "blue", lty = 2, lwd = 1)
abline(v = c(-1,1), col = "blue", lty = 2, lwd = 1)

# Volcano plot EdgeR
with(edge, plot(logFC, -log10(PValue),pch=1, main="Volcano plot: EdgeR",xlab="Log2FoldChange", ylab="Log10(pvalue)"))
#with(subset(edge, FDR>0.5), points(logFC, -log10(PValue), pch=20, col="gray"))
with(subset(edge, logFC < -1 & FDR<0.5), points(logFC, -log10(PValue), pch=1, col="red"))
with(subset(edge, logFC > 1 & FDR<0.5), points(logFC, -log10(PValue), pch=1, col="green"))
abline(h = 1, col = "blue", lty = 2, lwd = 1)
abline(v = c(-1,1), col = "blue", lty = 2, lwd = 1)


# Volcano plot ALDEx2
with(aldex, plot(diff.btw, -log10(we.ep), pch=1, main="Volcano plot: ALDEx2",xlab="Log2FoldChange", ylab="Log10(pvalue)"))
#with(subset(aldex, we.eBH>0.5), points(diff.btw, -log10(we.ep), pch=20, col="gray"))
with(subset(aldex, diff.btw < -1 & we.eBH<0.5), points(diff.btw, -log10(we.ep), pch=1, col="red"))
with(subset(aldex, diff.btw > 1 & we.eBH<0.5), points(diff.btw, -log10(we.ep), pch=1, col="green"))
abline(h = 1, col = "blue", lty = 2, lwd = 1)
abline(v = c(-1,1), col = "blue", lty = 2, lwd = 1)

dev.off()

# tables of subsets
deseq_up <- subset(deseq, log2FoldChange > 1 & padj<0.5)
edge_up <- subset(edge, logFC > 1 & FDR<0.5)
aldex_up <- subset(aldex, diff.btw < 1 & we.eBH<0.5)

X <- deseq$X
up <- merge(X, deseq_up$log2FoldChange, by="X")
up <- cbind(deseq_up$log2FoldChange, edge$logFC, aldex_up$diff.btw)

# to get the number of upregulated genes

deseq_up %>%
  ggplot(aes(x=X, y=log2FoldChange)) +
  geom_boxplot() +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("A boxplot with jitter") +
  xlab("")

# ------------------- STATISTICAL ANALYSIS ----------------------------




