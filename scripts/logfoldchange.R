# libraries
library("tidyverse")

# load the 3 log-fold change
deseq <- read.csv(unz("results/bulk/deseq2_DE.csv.zip", "deseq2_DE.csv"))
edge <- read.csv(unz("results/bulk/edgeR_DE.csv.zip", "edgeR_DE.csv"))
aldex <- read.csv("results/bulk/aldex2_DE.csv")



# Volcano plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=1, main="Volcano plot",xlab="Log2FoldChange", ylab="Log10(pvalue)"))
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=1, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=1, col="red"))
abline(h=-log10(.05), col="blue",lty=2)
abline(v=1, col="red",lty=2)
abline(v=-1, col="red",lty=2)
legend()
