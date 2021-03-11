library("tidyverse")
library("DESeq2")
library("SummarizedExperiment")

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

# Get expression matrix
expmatrix_dds <- vst(dds, fitType="local")
expmatrix <- assay(expmatrix_dds)


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

# Get top 30 high expressed genes
top <- order(rowMeans(expmatrix),decreasing=TRUE)[1:30]

# Heatmap
heatmap_data <- data.frame(sample_type = SummarizedExperiment::colData(dds)[,c("sample_type")], row.names = rownames(SummarizedExperiment::colData(dds)))
pheatmap::pheatmap(expmatrix[top,], cluster_rows=TRUE, show_rownames=TRUE, cluster_cols=TRUE, annotation_col=heatmap_data)

