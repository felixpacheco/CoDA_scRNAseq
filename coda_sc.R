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

ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
genes_names <- sc_counts$external_gene_name
IDs <- getBM(attributes=c('ensembl_gene_id',
                          'external_gene_name'),
             filters = 'mgi_symbol',
             values = genes_names,
             mart = ensembl)

# get annotated genes
annotated_genes_counts <- merge(IDs,sc_counts,by="external_gene_name")
annotated_genes_counts <- data.frame(column_to_rownames(annotated_genes_counts, var = "ensembl_gene_id"))
annotated_genes_counts <- annotated_genes_counts[,-1] 
genes <- row.names(annotated_genes_counts)

# ----------------- SingleCellExperiment -------------------------------------
# Now we need to read the SingleCell experiment
library("SingleCellExperiment")
require("scran")

### Make single cell experiment
sce <- SingleCellExperiment(assays = list(counts = annotated_genes_counts), 
                            colData = metadata)

# Pre doublet detection
#--- gene-annotation ---#
library(scater)
#rownames(sce) <- uniquifyFeatureNames(
 # rowData(sce)$gene, rowData(sce)$Symbol)

library(AnnotationHub)
library(org.Hs.eg.db)
library(AnnotationDbi)
ens.mm.v93 <- AnnotationHub()[["AH64461"]]


rowData(sce)$ensembl <- genes
rowData(sce)$SEQNAME <- mapIds(ens.mm.v93, keys=rowData(sce)$ensembl,
                                   keytype="GENEID", column="SEQNAME")

#--- quality-control ---#
is.mito <- rowData(sce)$SEQNAME == "MT"
stats <- perCellQCMetrics(sce, subsets=list(Mito=which(is.mito)))
qc <- quickPerCellQC(stats, percent_subsets="subsets_Mito_percent")
sce <- sce[,!qc$discard]

#--- normalization ---#
library(scran)
set.seed(101000110)
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters=clusters)
sce <- logNormCounts(sce)

#--- variance-modelling ---#
set.seed(00010101)
dec <- modelGeneVarByPoisson(sce)
top <- getTopHVGs(dec, prop=0.1)

#--- dimensionality-reduction ---#
library(BiocSingular)
set.seed(101010011)
sce <- denoisePCA(sce, technical=dec, subset.row=top)
sce <- runTSNE(sce, dimred="PCA")

#--- clustering ---#
snn.gr <- buildSNNGraph(sce, use.dimred="PCA", k=25)
colLabels(sce) <- factor(igraph::cluster_walktrap(snn.gr)$membership)

# Doublet detection
library("scDblFinder")
dbl.out <- findDoubletClusters(sce)
dbl.out


# PLOT , cannot plot there are not outliers
# library(scater)
# chosen.doublet <- rownames(dbl.out)[isOutlier(dbl.out$num.de, type="lower", log=TRUE)]
# library(scran)
# markers <- findMarkers(sce, direction="up")
# dbl.markers <- markers[[chosen.doublet]]
# 
# library(scater)
# chosen <- rownames(dbl.markers)[dbl.markers$Top <= 10]
# plotHeatmap(sce, order_columns_by="label", features=chosen, 
#             center=TRUE, symmetric=TRUE, zlim=c(-5, 5))
# ----------------------------------------------------------------------------
# ----------------- ALDEx2 package -------------------------------------------
# ----------------------------------------------------------------------------
# ----------------- ALDEx2 package -------------------------------------------
# ----------------------------------------------------------------------------
# the BioCManager packages need to be loaded always after the data wrangling because they interfere
library("ALDEx2")
conds <- metadata_counts$sample_type
x.all <- aldex(bulk_counts_clean, conds, mc.samples=128, test="t", effect=TRUE,include.sample.summary=FALSE, denom="all", verbose=FALSE)

# Write to file
write.csv(x.all, "aldex2_DE.csv", row.names = TRUE)

# clr
x <- aldex.clr(bulk_counts_clean, conds, mc.samples=128, denom="all", verbose=F)

x_df <- data.frame(x@analysisData)

# now reduce the columns to one
library(stringr)
nm1 <- str_replace(names(x_df),"\\..*", "")
x.new <- data.frame
x_df[paste0(unique(nm1), "_mean")] <- sapply(split.default(x_df, nm1), rowMeans)
x.new <- x_df[grep("_mean", names(x_df))] 
write.csv(x.new, "aldex2_clr.csv", row.names = TRUE)

# GLM
# one hot encode conds
# for(unique_value in unique(metadata_counts$sample_type)){
#   metadata_counts[paste("sample", unique_value, sep = ".")] <- ifelse(metadata_counts$sample_type == unique_value, 1, 0)
#   
# }
# o_conds <- metadata_counts[,c(9,10)]
# colnames(o_conds) <- c("T", "N")
# mm <- model.matrix(~ T + N, o_conds )
# clr <- aldex.clr(bulk_counts_clean, mm, mc.samples=128, denom="all", verbose=F)
# 
# glm.test <- aldex.glm(clr)
# glm.efftect <- aldex.glm.effect(clr, include.sample.summary=TRUE)
# 

#kw <- aldex.kw(x)
#kw.effect <- aldex.effect(x, include.sample.summary=TRUE, glm.conds=TRUE, CI=TRUE)
#print("kw and kw.effect done")
#write.csv(kw, gzfile("aldex2_kw.csv.gz"), row.names = TRUE)
#write.csv(kw.effect, gzfile("aldex2_kw_clr.csv.gz"), row.names = TRUE)
#print("files saved")





