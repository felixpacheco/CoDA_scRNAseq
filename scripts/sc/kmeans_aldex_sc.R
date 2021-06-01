rm(list=ls())
set.seed(123)
# Load the libraries
library("tidyverse")
library("xtable")
library(janitor)
library(cluster)
library(factoextra)
library(data.table)
library(stringr)
library(ggpubr)

metadata <- fread("/home/people/felpas/CoDA_scRNAseq/data/singlecell/metadata_sc.csv.gz", header=TRUE)
aldex_counts <- read.csv("/home/people/felpas/CoDA_scRNAseq/results/sc/aldex2_clr.csv.gz", sep = ",", header=TRUE)
print(head(aldex_counts))
### METADATA
# Get labels for plotting
metadata$rows <- rownames(metadata)
metadata$celltype <-  gsub("([A-Za-z]+).*", "\\1", metadata$cell_ontology_class)
metadata$label <- paste(metadata$celltype, metadata$rows, sep="_")

colnames(aldex_counts) <- substring(colnames(aldex_counts), 1,-5)
print(colnames(aldex_counts)[1])
# Replace colnames for metadata$labels for ALDEx2 values
labels_list_aldex <- c("gene_name")
columnnames_aldex <- colnames(aldex_counts)
for (i in columnnames_aldex) {
  new_element <- metadata[ metadata$cell == i, ]$label
  labels_list_aldex <- c(labels_list_aldex, new_element)
}
colnames(aldex_counts) <- labels_list_aldex

aldex_counts_t = as.data.frame(t(as.matrix(aldex_counts)))
aldex_counts_t <- clean_names(row_to_names(aldex_counts_t, row_number = 1))

print(head(aldex_counts_t))

### ---------- ALDEx2 DATA -----------
# Compute k-means
res_means_aldex = kmeans(aldex_counts_t, centers = 2, nstart = 20)

# Process output
res_viz_aldex <- as.data.frame(res_means_aldex$cluster)
names(res_viz_aldex)[1] <- "cluster"
res_viz_aldex <- cbind(sample_name = rownames(res_viz_aldex), res_viz_aldex)
rownames(res_viz_aldex) <- 1:nrow(res_viz_aldex)

res_viz_aldex$phenotype <- str_sub(res_viz_aldex$sample_name, 5)
res_viz_aldex$phenotype <- str_sub(res_viz_aldex$phenotype, start=1,end=1)

a1T <- nrow(res_viz_aldex[res_viz_aldex$cluster==1 & res_viz_aldex$phenotype=="f",])/nrow(res_viz_aldex[res_viz_aldex$cluster==1,]) 
a1N <- nrow(res_viz_aldex[res_viz_aldex$cluster==1 & res_viz_aldex$phenotype=="m",])/nrow(res_viz_aldex[res_viz_aldex$cluster==1,]) 
a2N <- nrow(res_viz_aldex[res_viz_aldex$cluster==2 & res_viz_aldex$phenotype=="m",])/nrow(res_viz_aldex[res_viz_aldex$cluster==2,]) 
a2T <- nrow(res_viz_aldex[res_viz_aldex$cluster==2 & res_viz_aldex$phenotype=="f",])/nrow(res_viz_aldex[res_viz_aldex$cluster==2,]) 

if (a1T>a1N) {
cat("ALDEx2 data -- Cluster purity:\n\tCluster 1:")
cat(a1T)
cat("\n\tCluster 2:")
cat(a2N)
} else {
cat("ALDEx2 data -- Cluster purity:\n\tCluster 1:")
cat(a1N)
cat("\n\tCluster 2:")
cat(a2T)
}


##---------- aldex -----------

aldex_countss <- as.data.frame(lapply(aldex_counts_t, as.numeric))
res.pca_aldex <- prcomp(aldex_countss, scale = FALSE)

df_pca_aldex <- as.data.frame(res.pca_aldex[["x"]])
labels <- row.names(aldex_counts_t)
df_pca_aldex$labels <- labels
df_pca_aldex$phenotype <- str_sub(df_pca_aldex$labels, start = 1, end = 1)

# Dimension reduction using PCA
res.pca_aldex <- prcomp(aldex_countss, scale = FALSE)
# Coordinates of individuals
ind.coord <- as.data.frame(get_pca_ind(res.pca_aldex)$coord)
# Add clusters obtained using the K-means algorithm
ind.coord$cluster <- factor(res_means_aldex$cluster)
# Add Species groups from the original data sett
ind.coord$label <- df_pca_aldex$phenotype
# Data inspection
head(ind.coord)

# Percentage of variance explained by dimensions
eigenvalue <- round(get_eigenvalue(res.pca_aldex), 1)
variance.percent <- eigenvalue$variance.percent
head(eigenvalue)
png(file="//home/people/felpas/CoDA_scRNAseq/results/sc/k_means_aldex1.png", height = 800, width = 1000, res=1000)
ggscatter(
  ind.coord, x = "Dim.1", y = "Dim.2", title = "ALDEx2 data clustering (k=2)",
  color = "cluster", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "label", size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", variance.percent[1], "% )" ),
  ylab = paste0("Dim 2 (", variance.percent[2], "% )" )
) +
  stat_mean(aes(color = cluster), size = 4)
dev.off()

png(file="//home/people/felpas/CoDA_scRNAseq/results/sc/pca_sc_aldex2.png", height = 800, width = 1000, res=1000)
ggplot(ind.coord, aes(x=Dim.1, y=Dim.2, shape=label, color=label)) +
  geom_point() + ggtitle("Principal component analysis projection") +
  xlab("Principal component 1") + ylab("Principal Compoment 2")
dev.off()
