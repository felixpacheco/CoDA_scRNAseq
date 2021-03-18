# -------------------- Load the libraries------------------------------------
library("tidyverse")
library("edgeR")
library("SummarizedExperiment")

# -------------------- Read the data ----------------------------------------
raw_bulk_counts <- read.csv("data/bulk/counts_bulk.tsv", sep = "\t")
metadata <- read.csv("data/bulk/metadata_bulk.tsv", sep = "\t")

# ---------------- Wrangle bulk counts --------------------------------------
bulk_counts <- raw_bulk_counts %>%
  subset(select = -Gene.Name) %>%
  column_to_rownames(var = "Gene.ID")

# -------------------- Clean metadata ---------------------------------------
# remove columns with all NA
metadata <- metadata %>% select_if(~ !all(is.na(.)))
# Rename columns, remove "Sample.Characteristic" and "." at the end of column name
names(metadata) <- gsub("Sample.Characteristic.", "", names(metadata))
names(metadata) <- sub(".$", "", names(metadata))
# Drop unused columns
# Age into only number
metadata <- metadata %>%
  rename(c("Ru" = "run", "sampling.site" = "sample_type")) %>%
  select(-contains("Ontology.Term")) %>%
  separate(age, c("age", "drop")) %>%
  subset(select = -drop)

# ----------------- EdgeR package -------------------------------------------
# Loading data to EdgeR
# Create EdgeR object
dgeFull <- DGEList(bulk_counts, group = metadata$sample_type)

# Add sample information to DGE object
dgeFull$metadata <- metadata
dgeFull

# Preparation for differential expression analysis
dge <- DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum) != 0, ],
               group=dgeFull$metadata$sample_type)

# Add sample information to DGE object
dge$metadata <- dgeFull$metadata

# DE analysis
# Calculate normalization factors + common&tagwise dispersions
dge <- calcNormFactors(dge)
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)
dge

# Exact test for the difference in expression between two sample types
dgeTest <- exactTest(dge)
dgeTest

# Obtain the FDR 
res <- topTags(dgeTest, n=nrow(dgeTest$table))
head(res)
res_df <- as.data.frame(res)

# ------------------ Data visualization --------------------------------------
# ------------------ Data dispersion -----------------------------------------

# ------------------- Volcano plot -------------------------------------------
# maybe use the library enhanced volcano to put some name tags to genes
# reset par and plot
par(mar = c(0, 0, 0, 0))
plot.new()
# Start our plot
par(mar = c(3, 3, 3, 3))
par(mfrow = c(1, 1))
# Make a basic volcano plot
with(res_df, plot(logFC, -log10(PValue), pch = 1, main = "Volcano plot", xlab = "Log2FoldChange", ylab = "Log10(pvalue)"))
with(subset(res_df, FDR < .05), points(logFC, -log10(PValue), pch = 1, col = "blue"))
with(subset(res_df, FDR < .05 & abs(logFC) > 1), points(logFC, -log10(PValue), pch = 1, col = "red"))
abline(h = -log10(.05), col = "blue", lty = 2)
abline(v = 1, col = "red", lty = 2)
abline(v = -1, col = "red", lty = 2)
# legend()

# ------------------- MA Plot ------------------------------------------------
# DESeq intregrated function
plotMA(res)
