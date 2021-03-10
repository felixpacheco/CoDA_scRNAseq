library("tidyverse")
library("DESeq2")

# Read the data
raw_bulk_counts <- read.csv("data/bulk/counts_bulk.tsv", sep = "\t")
metadata <- read.csv("data/bulk/metadata_bulk.tsv", sep = "\t")

# Wrangle bulk counts
bulk_counts <- raw_bulk_counts %>% subset(select=-Gene.Name)

# Clean metadata
# remove columns with all NA
metadata <- metadata %>% select_if(~ !all(is.na(.)))
# Rename columns, remove "Sample.Characteristic" and "." at the end of column name
names(metadata) <- gsub("Sample.Characteristic.", "", names(metadata))
names(metadata) <- sub(".$", "", names(metadata))
# Drop unused columns
# Age into only number
metadata <- metadata %>%
  rename(c("Ru" = "run","sampling.site"="sample_type")) %>%
  select(-contains("Ontology.Term")) %>%
  separate(age, c("age", "drop")) %>%
  subset(select = -drop)

# Loading data to DESeq2
# Create DESeq object
dds <- DESeqDataSetFromMatrix(countData=bulk_counts, 
                              colData=metadata, 
                              design=~sample_type, tidy = TRUE)

dds <- DESeq(dds)
res <- results(dds)
head(results(dds, tidy=TRUE)) 

summary(res) #summary of results

res_ord <- res[order(res$padj),]
head(res_ord)


par(mfrow=c(2,3))
plotCounts(dds, gene="ENSG00000133636", intgroup="sample_type")
plotCounts(dds, gene="ENSG00000000005", intgroup="sample_type")
plotCounts(dds, gene="ENSG00000000938", intgroup="sample_type")
plotCounts(dds, gene="ENSG00000000460", intgroup="sample_type")
plotCounts(dds, gene="ENSG00000132872", intgroup="sample_type")
plotCounts(dds, gene="ENSG00000096088", intgroup="sample_type")

# Volcano plot
#reset par
par(mfrow=c(1,1))
# Make a basic volcano plot
with(res, plot(log2FoldChange, -log10(pvalue), pch=1, main="Volcano plot",xlab="Log2FoldChange", ylab="Log10(pvalue)"))
# Add colored points: blue if padj<0.01, red if log2FC>1 and padj<0.05)
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=1, col="blue"))
with(subset(res, padj<.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=1, col="red"))
abline(h=-log10(.05), col="blue",lty=2)
abline(h=log10(.05), col="blue",lty=2)
abline(v=2, col="red",lty=2)
abline(v=-2, col="red",lty=2)



