# -------------------- Load the libraries------------------------------------
library("tidyverse")
library("xtable")

# -------------------- Read the data ----------------------------------------
sc_counts <- read_csv("data/singlecell/counts_sc.csv.gz")
metadata <- read_csv("data/singlecell/metadata_sc.csv.gz")

# Clean and wrangle
metadata <- metadata[,-1]
names(sc_counts)[names(sc_counts) == 'Row.names'] <- 'gene'
sc_counts <- data.frame(column_to_rownames(sc_counts, var = "gene"))
sc_counts <- sc_counts[,-1] 
# take the rows with all 0s 
sc_counts <- sc_counts[apply(sc_counts[, -1], 1, function(x) !all(x == 0)), ]

# ----------------- SingleCellExperiment -------------------------------------
# Now we need to read the SingleCell experiment
library("SingleCellExperiment")
require("scran")

### Make single cell experiment
sce <- SingleCellExperiment(assays = list(counts = sc_counts), 
                            colData = metadata)

# Pre doublet detection
#--- gene-annotation ---#
library(scater)
rownames(sce) <- uniquifyFeatureNames(
  rowData(sce)$Ensembl, rowData(sce)$Symbol)

library(AnnotationHub)
ens.mm.v93 <- AnnotationHub()[["AH64461"]]
library(data.table)
df$HeaderName <- row.names()
rowData(sce)$SEQNAME <- mapIds(ens.mm.v93, keys=rowData(sce)$Ensembl,
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
sce.mam <- denoisePCA(sce, technical=dec, subset.row=top)
sce.mam <- runTSNE(sce, dimred="PCA")

#--- clustering ---#
snn.gr <- buildSNNGraph(sce.mam, use.dimred="PCA", k=25)
colLabels(sce) <- factor(igraph::cluster_walktrap(snn.gr)$membership)

# Doublet detection
library("scDblFinder")
dbl.out <- findDoubletClusters(sce)
dbl.out


# ----------------------------------------------------------------------------
# ----------------- ALDEx2 package -------------------------------------------
# ----------------------------------------------------------------------------



