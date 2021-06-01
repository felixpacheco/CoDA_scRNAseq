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

#### ----- READ THE DATA -----
raw_counts <- read.csv("/home/people/felpas/CoDA_scRNAseq/data/singlecell/counts_sc.csv.gz", sep = ",")
metadata <- fread("/home/people/felpas/CoDA_scRNAseq/data/singlecell/metadata_sc.csv.gz", header=TRUE)

clr_counts <- read.csv("/home/people/felpas/clr_sc.csv.gz", sep = ",")
deseq_counts <- read.csv("/home/people/felpas/SC_deseq2_DE.csv.gz", sep = ",", header=TRUE)
edge_counts <- read.csv("/home/people/felpas/SC_edgeR_DE.csv.gz", sep = ",", header=TRUE)
print("loaded data succesfully")

### ----- WRANGLE RAW DATA -----
### METADATA
# Get labels for plotting
metadata$rows <- rownames(metadata)
metadata$celltype <-  gsub("([A-Za-z]+).*", "\\1", metadata$cell_ontology_class)
metadata$label <- paste(metadata$celltype, metadata$rows, sep="_")

### RAW COUNTS
# Drop X column and remove genes with 0 counts everywhere
raw_counts = subset(raw_counts, select = -c(X))
raw_counts<- raw_counts[apply(raw_counts[,-1], 1, function(x) !all(x==0)),]

# Replace colnames for metadata$labels
labels_list <- c("gene_name")
columnnames_raw <- colnames(raw_counts)
for (i in columnnames_raw) {
  new_element <- metadata[ metadata$cell == i, ]$label
  labels_list <- c(labels_list, new_element)
}
colnames(raw_counts) <- labels_list

### DESEQ AND EDGE : Delete non relevant columns
deseq_counts <- subset(deseq_counts, select = -c(baseMean,log2FoldChange,lfcSE, stat,pvalue,padj))
edge_counts <- subset(edge_counts, select = -c(logCPM,PValue,dispersion, logFC,LR,FDR))

# Replace colnames for metadata$labels for deseq2
labels_list_deseq <- c("gene_name")
columnnames_deseq <- colnames(deseq_counts)
for (i in columnnames_deseq) {
  new_element <- metadata[ metadata$cell == i, ]$label
  print(new_element)
  labels_list_deseq <- c(labels_list_deseq, new_element)
}
colnames(deseq_counts) <- labels_list_deseq

# Replace colnames for metadata$labels for edge R
labels_list_edge <- c("gene_name")
columnnames_edge <- colnames(edge_counts)
for (i in columnnames_edge) {
  new_element <- metadata[ metadata$cell == i, ]$label
  labels_list_edge <- c(labels_list_edge, new_element)
}
colnames(edge_counts) <- labels_list_edge

# Replace colnames for metadata$labels for CLR values
clr_counts <- subset(clr_counts, select=-c(Unnamed..0))
labels_list_clr <- c("gene_name")
columnnames_clr <- colnames(clr_counts)
for (i in columnnames_clr) {
  new_element <- metadata[ metadata$cell == i, ]$label
  labels_list_clr <- c(labels_list_clr, new_element)
}
colnames(clr_counts) <- labels_list_clr

# Transpose counts matrix

raw_counts_t = as.data.frame(t(as.matrix(raw_counts)))
raw_counts_t <- clean_names(row_to_names(raw_counts_t, row_number = 1))

# Transpose counts matrix

deseq_counts_t = as.data.frame(t(as.matrix(deseq_counts)))
deseq_counts_t <- clean_names(row_to_names(deseq_counts_t, row_number = 1))

# Transpose counts matrix

edge_counts_t = as.data.frame(t(as.matrix(edge_counts)))
edge_counts_t <- clean_names(row_to_names(edge_counts_t, row_number = 1))

# Transpose counts matrix

clr_counts_t = as.data.frame(t(as.matrix(clr_counts)))
clr_counts_t <- clean_names(row_to_names(clr_counts_t, row_number = 1))

print("wrangled data succesfully")

## KNN algorithm

# RAW DATA
set.seed(1996)
ran <- sample(1:nrow(raw_counts_t), 0.9 * nrow(raw_counts_t)) 
 
##extract training set
#bulk_counts <- bulk_counts_t[,!names(bulk_counts_t) %in% "Ensembl_gene_id"]
raw_counts_train <- raw_counts_t[ran,] 

##extract testing set
raw_counts_test <- raw_counts_t[-ran,] 

## extract labels
targets_train<- str_sub(rownames(raw_counts_t), start=1,end=1)[ran]
targets_test <-  str_sub(rownames(raw_counts_t), start=1,end=1)[-ran]
# Compute KNN function
library(class)

accuracy <- function(x){sum(diag(x)/(sum(rowSums(x)))) * 100}

k_to_try = 1:15
err_k_raw = rep(x = 0, times = length(k_to_try))

for (i in seq_along(k_to_try)) {
  pred = knn(train = raw_counts_train, 
             test  = raw_counts_test, 
             cl    = targets_train, 
             k     = k_to_try[i])
  tab = table(pred, targets_test)
  err_k_raw[i] <- accuracy(tab)
}
print("KNN raw computed")

# DESeq2 DATA
ran <- sample(1:nrow(deseq_counts_t), 0.9 * nrow(deseq_counts_t)) 
 
##extract training set
#deseq_counts <- deseq_counts[,!names(deseq_counts) %in% "X"]
deseq_counts_train <- deseq_counts_t[ran,] 

##extract testing set
deseq_counts_test <- deseq_counts_t[-ran,] 

## extract labels
targets_train<- str_sub(rownames(deseq_counts_t), start=1,end=1)[ran]
targets_test <-  str_sub(rownames(deseq_counts_t), start=1,end=1)[-ran]
# Compute KNN function
library(class)

accuracy <- function(x){sum(diag(x)/(sum(rowSums(x)))) * 100}

k_to_try = 1:15
err_k_deseq = rep(x = 0, times = length(k_to_try))

for (i in seq_along(k_to_try)) {
  pred = knn(train = deseq_counts_train, 
             test  = deseq_counts_test, 
             cl    = targets_train, 
             k     = k_to_try[i])
  tab = table(pred, targets_test)
  err_k_deseq[i] <- accuracy(tab)
}
print("KNN deseq computed")
# EdgeR
ran <- sample(1:nrow(edge_counts_t), 0.9 * nrow(edge_counts_t)) 
 
##extract training set
edge_counts <- edge_counts_t[,!names(edge_counts_t) %in% "X"]
edge_counts_train <- edge_counts_t[ran,] 

##extract testing set
edge_counts_test <- edge_counts_t[-ran,] 

## extract labels
targets_train<- str_sub(rownames(edge_counts_t), start=1,end=1)[ran]
targets_test <-  str_sub(rownames(edge_counts_t), start=1,end=1)[-ran]
# Compute KNN function
library(class)

accuracy <- function(x){sum(diag(x)/(sum(rowSums(x)))) * 100}

k_to_try = 1:15
err_k_edge = rep(x = 0, times = length(k_to_try))

for (i in seq_along(k_to_try)) {
  pred = knn(train = edge_counts_train, 
             test  = edge_counts_test, 
             cl    = targets_train, 
             k     = k_to_try[i])
  tab = table(pred, targets_test)
  err_k_edge[i] <- accuracy(tab)
}
print("KNN edge computed")
# CLR
ran <- sample(1:nrow(clr_counts_t), 0.9 * nrow(clr_counts_t)) 
 
##extract training set
clr_counts <- clr_counts_t[,!names(clr_counts_t) %in% "X"]
clr_counts_train <- clr_counts_t[ran,] 

##extract testing set
clr_counts_test <- clr_counts_t[-ran,] 

## extract labels
targets_train<- str_sub(rownames(clr_counts_t), start=1,end=1)[ran]
targets_test <-  str_sub(rownames(clr_counts_t), start=1,end=1)[-ran]
# Compute KNN function
library(class)

accuracy <- function(x){sum(diag(x)/(sum(rowSums(x)))) * 100}

k_to_try = 1:15
err_k_clr = rep(x = 0, times = length(k_to_try))

for (i in seq_along(k_to_try)) {
  pred = knn(train = clr_counts_train, 
             test  = clr_counts_test, 
             cl    = targets_train, 
             k     = k_to_try[i])
  tab = table(pred, targets_test)
  err_k_clr[i] <- accuracy(tab)
}
print("KNN clr computed")

png(file="/home/people/felpas/CoDA_scRNAseq/results/sc/raw_KNN_sc_15.png", height = 800, width = 1000, res=1000)
plot(err_k_raw, type = "b", col = "darkolivegreen4",
     lty = 1, xlab = "k, number of neighbors", ylab = "Accuracy",
     main = "Accuracy with k-neighbors with raw data",ylim=c(60,103))
dev.off()

png(file="/home/people/felpas/CoDA_scRNAseq/results/sc/deseq_KNN_sc_15.png", height = 800, width = 1000, res=1000)
plot(err_k_deseq, type = "b", col = "darkolivegreen4",
     lty = 1, xlab = "k, number of neighbors", ylab = "Accuracy",
     main = "Accuracy with k-neighbors with DESeq2",ylim=c(60,103))
dev.off()

png(file="/home/people/felpas/CoDA_scRNAseq/results/sc/edge_KNN_sc_15.png", height = 800, width = 1000, res=1000)
plot(err_k_edge, type = "b", col = "darkolivegreen4",
     lty = 1, xlab = "k, number of neighbors", ylab = "Accuracy",
     main = "Accuracy with k-neighbors with Edge R",ylim=c(20,103))
dev.off()

png(file="//home/people/felpas/CoDA_scRNAseq/results/sc/clr_KNN_sc_15.png", height = 800, width = 1000, res=1000)
plot(err_k_clr, type = "b", col = "darkolivegreen4",
     lty = 1, xlab = "k, number of neighbors", ylab = "Accuracy",
     main = "Accuracy with k-neighbors with CLR",ylim=c(60,103))
dev.off()

print("graphs computed normally")