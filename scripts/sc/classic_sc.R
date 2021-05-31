rm(list=ls())
options(warn=-1)
# -------------------- Load the libraries------------------------------------
library("tidyverse")
library("xtable")

# -------------------- Read the data ----------------------------------------
sc_counts <- read_csv("data/singlecell/counts_sc.csv.gz")
metadata <- read_csv("data/singlecell/metadata_sc.csv.gz")

# Clean and wrangle
metadata <- metadata[,-1]
names(sc_counts)[names(sc_counts) == 'Row.names'] <- 'external_gene_name'

sc_counts <- sc_counts[,-1] 
# take the rows with all 0s 
sc_counts <- sc_counts[apply(sc_counts[, -1], 1, function(x) !all(x == 0)), ]

# get ensmbl gene ids
library(biomaRt)

#ensembl <- useMart("ensembl", host="uswest.ensembl.org", dataset="mmusculus_gene_ensembl", ensemblRedirect = FALSE)
genes_names <- sc_counts$external_gene_name

#IDs <- getBM(attributes=c('ensembl_gene_id',
#                          'external_gene_name'),
#             filters = 'mgi_symbol',
#             values = genes_names,
#             mart = ensembl)
IDs <- read_csv("data/martquery_0530110047_39.txt")
names(IDs)[names(IDs) == 'Gene stable ID'] <- 'ensembl_gene_id'
names(IDs)[names(IDs) == 'Gene name'] <- 'external_gene_name'


# get annotated genes
annotated_genes_counts <- merge(IDs,sc_counts,by="external_gene_name")
annotated_genes_counts <- data.frame(column_to_rownames(annotated_genes_counts, var = "ensembl_gene_id"))
annotated_genes_counts <- annotated_genes_counts[,-1] 
genes <- row.names(annotated_genes_counts)

# ----------------- SingleCellExperiment -------------------------------------
# Now we need to read the SingleCell experiment
require("SingleCellExperiment")
require("scran")
# require(AnnotationHub)
# require(org.Hs.eg.db)
# require(AnnotationDbi)
# require(scater)
# require(BiocSingular)
# require("scDblFinder")
# 
### Make single cell experiment
sce <- SingleCellExperiment(assays = list(counts = annotated_genes_counts), 
                            colData = metadata)
# 
# # Pre doublet detection
# #--- gene-annotation ---#
# #rownames(sce) <- uniquifyFeatureNames(
# # rowData(sce)$gene, rowData(sce)$Symbol)
# 
# ens.mm.v93 <- AnnotationHub()[["AH64461"]]
# 
# rowData(sce)$ensembl <- genes
# rowData(sce)$SEQNAME <- mapIds(ens.mm.v93, keys=rowData(sce)$ensembl,
#                                keytype="GENEID", column="SEQNAME")
# 
# #--- quality-control ---#
# is.mito <- rowData(sce)$SEQNAME == "MT"
# stats <- perCellQCMetrics(sce, subsets=list(Mito=which(is.mito)))
# qc <- quickPerCellQC(stats, percent_subsets="subsets_Mito_percent")
# sce <- sce[,!qc$discard]
# 
# #--- normalization ---#
# set.seed(101000110)
# clusters <- quickCluster(sce)
# sce <- computeSumFactors(sce, clusters=clusters)
# sce <- logNormCounts(sce)
# 
# #--- variance-modelling ---#
# set.seed(00010101)
# dec <- modelGeneVarByPoisson(sce)
# top <- getTopHVGs(dec, prop=0.1)
# 
# #--- dimensionality-reduction ---#
# set.seed(101010011)
# sce <- denoisePCA(sce, technical=dec, subset.row=top)
# sce <- runTSNE(sce, dimred="PCA")
# 
# #--- clustering ---#
# snn.gr <- buildSNNGraph(sce, use.dimred="PCA", k=25)
# colLabels(sce) <- factor(igraph::cluster_walktrap(snn.gr)$membership)
# 
# # Doublet detection
# dbl.out <- findDoubletClusters(sce)
# dbl.out

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

dds_mat <-estimateSizeFactors(dds_mat, type = 'iterate')
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


# ----------------------------------------------------------------------------
# ----------------- ALDEx2 package -------------------------------------------
# ----------------------------------------------------------------------------
# no use of Single Cell Experiment
# the BioCManager packages need to be loaded always after the data wrangling because they interfere
library("ALDEx2")
conds <- metadata$cell_ontology_class
x.all <- aldex(annotated_genes_counts, conds, mc.samples=128, test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)

# Write to file
write.csv(x.all, "aldex2_DE.csv", row.names = TRUE)

# clr
x <- aldex.clr(annotated_genes_counts, conds, mc.samples=128, denom="all", verbose=F)

x_df <- data.frame(x@analysisData)

# now reduce the columns to one
library(stringr)
nm1 <- str_replace(names(x_df),"\\..*", "")
x.new <- data.frame
x_df[paste0(unique(nm1), "_mean")] <- sapply(split.default(x_df, nm1), rowMeans)
x.new <- x_df[grep("_mean", names(x_df))] 
write.csv(x.new, "aldex2_clr.csv", row.names = TRUE)

