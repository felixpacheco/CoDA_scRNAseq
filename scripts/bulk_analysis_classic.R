# -------------------- Load the libraries------------------------------------
library("tidyverse")
library("DESeq2")
library("edgeR")
library("SummarizedExperiment")
# -------------------- Read the data ----------------------------------------
raw_bulk_counts <- read.csv("data/bulk/counts_bulk.tsv", sep = "\t")
metadata <- read.csv("data/bulk/metadata_bulk.tsv", sep = "\t")
# ---------------- Wrangle bulk counts --------------------------------------
bulk_counts <- raw_bulk_counts %>% subset(select = -Gene.Name)
bulk_counts_edgeR <- raw_bulk_counts %>%
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
  rename(c("run" = "Ru", "sample_type" = "sampling.site")) %>%
  select(-contains("Ontology.Term")) %>%
  separate(age, c("age", "drop")) %>%
  subset(select = -drop)
# ----------------------------------------------------------------------------
# ----------------- DESeq2 package -------------------------------------------
# ----------------------------------------------------------------------------
# Loading data to DESeq2
# Create DESeq object
dds_mat <- DESeqDataSetFromMatrix(
  countData = bulk_counts,
  colData = metadata,
  design = ~sample_type, tidy = TRUE
)
dds <- DESeq(dds_mat)
res <- results(dds)

head(results(dds, tidy = TRUE))
summary(res) # summary of results

vsd <- vst(dds_mat, blind = FALSE)
res_vsd <- cbind(as.data.frame(res), assay(vsd))
dge_dds <- as.data.frame(res_vsd[order(res_vsd$pvalue), ])
head(dge_dds, 10)
# ----------------------------------------------------------------------------
# ----------------- EdgeR package -------------------------------------------
# ----------------------------------------------------------------------------
# Loading data to EdgeR
# Create EdgeR object
dgeFull <- DGEList(bulk_counts_edgeR, group = metadata$sample_type)

# Add sample information to DGE object
dgeFull$metadata <- metadata
dgeFull

# Preparation for differential expression analysis
dge <- DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum) != 0, ],
  group = dgeFull$metadata$sample_type
)
# Add sample information to DGE object
dge$metadata <- dgeFull$metadata

# DE analysis
# Calculate normalization factors + common&tagwise dispersions
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge)
# dge <- estimateCommonDisp(dge)
# dge <- estimateTagwiseDisp(dge)
dge

# Exact test for the difference in expression between two sample types
# dgeTest <- exactTest(dge)
# dgeTest
# Fit glm function
sample_type_mat <- relevel(factor(metadata$sample_type), ref = "tumor tissue")
edesign <- model.matrix(~sample_type_mat)
fit <- glmFit(dge, edesign)
lrt <- glmLRT(fit)
res_edgeR <- as.data.frame(topTags(lrt, n = Inf))
res_edgeR$dispersion <- lrt$dispersion
res_edgeR <- merge(res_edgeR, lrt$fitted.values, by = "row.names")
rownames(res_edgeR) <- res_edgeR$Row.names
res_edgeR$Row.names <- NULL
res_edgeR <- res_edgeR[order(res_edgeR$PValue), ]
head(res_edgeR, 10)

# ----------------- Outfiles -----------------------------
# Write to file
write.csv(dge_dds, "deseq2_DE.csv", row.names = TRUE)
# Write to file
write.csv(res_edgeR, "edgeR_DE.csv", row.names = TRUE)
