# -------------------- Load the libraries------------------------------------
library("tidyverse")
library("xtable")


# -------------------- Read the data ----------------------------------------
sc_counts <- read_csv("data/singlecell/counts_sc.csv.gz")
metadata <- read_csv("data/singlecell/metadata_sc.csv.gz")

# Clean and wrangle
metadata <- metadata[,-1]
sc_counts <- data.frame(column_to_rownames(sc_counts, var = "Row.names"))
sc_counts <- sc_counts[,-1] 

# ----------------------------------------------------------------------------
# ----------------- EdgeR package --------------------------------------------
# ----------------------------------------------------------------------------
# Now we need to read the SingleCell experiment
library("SingleCellExperiment")
library("edgeR")
require("scran")

### Make single cell experiment
sce <- SingleCellExperiment(assays = sc_counts, 
                            colData = metadata)

# Now proceed with EdgeR
# Loading data to EdgeR - Create EdgeR object
dgeFull <- convertTo(sce, type="edgeR", assay.type = 1)

# Preparation for differential expression analysis
dge <- DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum) != 0, ],
               group = dgeFull$metadata$cell_ontology_class
)
# Add sample information to DGE object
dge$metadata <- dgeFull$metadata

# DE analysis
# Calculate normalization factors + common&tagwise dispersions
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge)
dge

# Fit glm function for log-likelihood ratio test (LTR)
sample_type_mat <- relevel(factor(metadata$cell_ontology_class), ref = "tumor tissue")
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








# ----------------------------------------------------------------------------
# ----------------- DESeq2 package -------------------------------------------
# ----------------------------------------------------------------------------
library("DESeq2")
# Loading data to DESeq2
# Create DESeq object
dds_mat <- DESeqDataSetFromMatrix(
  countData = bulk_counts_clean,
  colData = metadata_counts,
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
# ----------------- EdgeR package --------------------------------------------
# ----------------------------------------------------------------------------
# Wrangle for edge R
row.names(bulk_counts_clean) <- bulk_counts_clean$Ensembl_gene_id
bulk_counts_clean <- bulk_counts_clean[,-1] 

library("edgeR")
# Loading data to EdgeR
# Create EdgeR object
dgeFull <- DGEList(bulk_counts_clean, group = metadata_counts$sample_type)

# Add sample information to DGE object
dgeFull$metadata <- metadata_counts
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
# Fit glm function for log-likelihood ratio test (LTR)
sample_type_mat <- relevel(factor(metadata_counts$sample_type), ref = "tumor tissue")
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
