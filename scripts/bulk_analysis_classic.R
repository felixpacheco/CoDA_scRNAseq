# -------------------- Load the libraries------------------------------------
library("tidyverse")
library("xtable")
# -------------------- Read the data ----------------------------------------
raw_bulk_counts <- read.csv("data/bulk/counts_bulk_patient.tsv", sep = "\t")
metadata_raw <- read.csv("data/bulk/metadata_bulk.tsv", sep = "\t")
# ---------------- Wrangle bulk counts --------------------------------------
bulk_counts_edge <- raw_bulk_counts %>%
  column_to_rownames(var = "Ensembl_gene_id")
bulk_counts <- raw_bulk_counts
bulk_counts[63152,"Ensembl_gene_id"] <-"TC_"
bulk_counts <- bulk_counts[-c(63152), ] 
bulk_counts_edge <- bulk_counts_edge[-c(63152), ] 
# Import head of raw bulk counts to latex table
xtable(bulk_counts_edge[c(1:5),c(1:4, 39)])

# -------------------- Clean metadata ---------------------------------------
# remove columns with all NA
metadata <- metadata_raw %>% select_if(~ !all(is.na(.)))
# Rename columns, remove "Sample.Characteristic" and "." at the end of column name
names(metadata) <- gsub("Sample.Characteristic.", "", names(metadata))
names(metadata) <- sub(".$", "", names(metadata))
# Drop unused columns
# Age into only number
metadata <- metadata %>%
  rename("sample_type" = "sampling.site") %>%
  select(-contains("Ontology.Term")) %>%
  separate(age, c("age", "drop")) %>%
  subset(select = -c(Ru,drop, organism, Factor.Value.sampling.site, Analyse))
# Now get only metadata per individual not per run
metadata <- metadata %>% mutate(sample = ifelse(sample_type=='tumor tissue', paste(individual,"T", sep=""), paste(individual,"N", sep="")))

metadata_counts <- metadata %>% group_by(sample) %>% 
  distinct
# ----------------------------------------------------------------------------
# ----------------- DESeq2 package -------------------------------------------
# ----------------------------------------------------------------------------
library("DESeq2")
# Loading data to DESeq2
# Create DESeq object
dds_mat <- DESeqDataSetFromMatrix(
  countData = bulk_counts,
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
library("edgeR")
# Loading data to EdgeR
# Create EdgeR object
dgeFull <- DGEList(bulk_counts_edge, group = metadata_counts$sample_type)

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
# Fit glm function
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
