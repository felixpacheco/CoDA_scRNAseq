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
# take the rows with all 0s 
sc_counts <- sc_counts[apply(sc_counts[, -1], 1, function(x) !all(x == 0)), ]

# ----------------------------------------------------------------------------
# ----------------- Single Cell Experiment package ---------------------------
# ----------------------------------------------------------------------------
# Now we need to read the SingleCell experiment
library("SingleCellExperiment")
require("scran")

### Make single cell experiment
sce <- SingleCellExperiment(assays = sc_counts, 
                            colData = metadata)

# ----------------------------------------------------------------------------
# ----------------- EdgeR package --------------------------------------------
# ----------------------------------------------------------------------------
library("edgeR")
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
sample_type_mat <- relevel(factor(metadata$cell_ontology_class), ref = "fibroblast")
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

# Now proceed with DESeq2
# Loading data to DESeq2 - Create DESeq2 object
dds_mat <- convertTo(sce, type="DESeq2", assay.type = 1)

# Create DESeq object
dds <- DESeq(dds_mat)
res <- results(dds)

head(results(dds, tidy = TRUE))
summary(res) # summary of results

vsd <- vst(dds_mat, blind = FALSE)
res_vsd <- cbind(as.data.frame(res), assay(vsd))
dge_dds <- as.data.frame(res_vsd[order(res_vsd$pvalue), ])
head(dge_dds, 10)

# ------------------------ Outfiles ------------------------------------------
# Write to file
write.csv(dge_dds,  gzfile("SC_deseq2_DE.csv.gz"), row.names = TRUE)
# Write to file
write.csv(res_edgeR, gzfile("SC_edgeR_DE.csv.gz"), row.names = TRUE)




library("tidyverse")
df_edge <- read_csv("SC_edgeR_DE.csv.gz")
df_deseq <- read_csv("SC_deseq2_DE.csv.gz")