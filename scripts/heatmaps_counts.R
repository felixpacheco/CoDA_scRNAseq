# -------------------- Load the libraries------------------------------------
library("tidyverse")
library("xtable")
library("viridis")
library("hrbrthemes")
# -------------------- Read the data ----------------------------------------
raw_bulk_counts <- read.csv("data/bulk/counts_bulk_patient.tsv", sep = "\t")
metadata_raw <- read.csv("data/bulk/metadata_bulk.tsv", sep = "\t")
# ---------------- Wrangle bulk counts --------------------------------------
bulk_counts <- raw_bulk_counts
bulk_counts[63152, "Ensembl_gene_id"] <- "TC_"
bulk_counts <- bulk_counts[-c(63152), ]
bulk_counts_clean <- bulk_counts[apply(bulk_counts[, -1], 1, function(x) !all(x == 0)), ]

# Import head of raw bulk counts to latex table
xtable(bulk_counts_clean[c(1:5), c(1:4, 39)])

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
  subset(select = -c(Ru, drop, organism, Factor.Value.sampling.site, Analyse))
# Now get only metadata per individual not per run
metadata <- metadata %>% mutate(sample = ifelse(sample_type == "tumor tissue", paste(individual, "T", sep = ""), paste(individual, "N", sep = "")))

metadata_counts <- metadata %>%
  group_by(sample) %>%
  distinct()

# ---------------------- start general visualitzations -------------------------
# Wrangle for visualizations
row.names(bulk_counts_clean) <- bulk_counts_clean$Ensembl_gene_id
bulk_counts_clean <- bulk_counts_clean[, -1]

# heat map raw counts random sample
set.seed(123539)
cols <- sample(1:ncol(bulk_counts_clean), 50)
rows <- sample(1:nrow(bulk_counts_clean), 50)

# First pick random cols and rows
bulk_counts_clean <- bulk_counts_clean.sort_index()
random_raw_counts <- bulk_counts_clean[rows, cols]

# We need a long dataframe to plot the heatmap
long_random_raw_counts <- random_raw_counts %>%
  rownames_to_column() %>%
  gather(colname, Counts, -rowname)
head(long_random_raw_counts)

long_random_raw_counts <- long_random_raw_counts %>% mutate(norm = (Counts - min(Counts)) / (max(Counts) - min(Counts)))
# Now create the normalized 0-1 column:

png(file = "heatmap_raw_counts.png", height = 7, width = 10, units = "in", res = 300)
p1 <- ggplot(long_random_raw_counts, aes(x = colname, y = rowname, fill = norm)) +
  geom_tile() +
  labs(x = "Patient ID", y = "Gene ID", tag = "A", fill = "Counts") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "none", text = element_text(size = 20)
  ) +
  scale_fill_viridis(discrete = FALSE) #+
  #scale_y_discrete(breaks=long_random_raw_counts$rowname[seq(1,length(long_random_raw_counts$rowname),by=2)]) +
  #scale_x_discrete(breaks=long_random_raw_counts$colname[seq(1,length(long_random_raw_counts$colname),by=10)])

dev.off()

# what if we do the same after DESEQ normalization and EDGER
# open the files
dds_norm <- read.csv(unz("results/bulk/deseq2_DE.csv.zip", "deseq2_DE.csv"))
edge_norm <- read.csv(unz("results/bulk/edgeR_DE.csv.zip", "edgeR_DE.csv"))

# get top 500 and write to file
# sorted_by_log_fc <- dds_norm %>% arrange(log2FoldChange)
# Write to file
# write.csv(head(dds_norm,500), "deseq2_top500.csv", row.names = TRUE)

# do the same plot as before:
# sort the same
dds_norm <- dds_norm %>% arrange(X)
edge_norm <- edge_norm %>% arrange(X)

# clean data
row.names(dds_norm) <- dds_norm$X
dds_norm <- dds_norm[, -1]
row.names(edge_norm) <- edge_norm$X
edge_norm <- edge_norm[, -1]

dds_norm_counts <- dds_norm[, -c(1, 2, 3, 4, 5, 6)]
edge_norm_counts <- edge_norm[, -c(1, 2, 3, 4, 5, 6)]

# Import head of raw bulk counts to latex table
xtable(dds_norm_counts[c(1:5), c(1:4, 39)])
xtable(edge_norm_counts[c(1:5), c(1:4, 39)])
# heat map normalized counts random sample
# get random cols and rows
random_dds_norm <- dds_norm_counts[rows, cols]
random_edge_norm <- edge_norm_counts[rows, cols]

# We need a long dataframe to plot the heatmap
long_random_edge_norm <- random_edge_norm %>%
  rownames_to_column() %>%
  gather(colname, Counts, -rowname)
head(long_random_edge_norm)

long_random_dds_norm <- random_dds_norm %>%
  rownames_to_column() %>%
  gather(colname, Counts, -rowname)
head(long_random_dds_norm)

# Now create the normalized 0-1 column:
long_random_dds_norm <- long_random_dds_norm %>% mutate(norm = (Counts - min(Counts)) / (max(Counts) - min(Counts)))
long_random_edge_norm <- long_random_edge_norm %>% mutate(norm = (Counts - min(Counts)) / (max(Counts) - min(Counts)))

png(file = "heatmap_dds_norm_counts.png", height = 7, width = 10, units = "in", res = 300)
p2 <- ggplot(long_random_dds_norm, aes(x = colname, y = rowname, fill = norm)) +
  geom_tile() +
  labs(x = "Patient ID", y = "Gene ID", tag = "B", fill = "Counts") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    text = element_text(size = 20)
  ) +
  scale_fill_viridis(discrete = FALSE) #+
  #scale_x_discrete(breaks=long_random_dds_norm$colname[seq(1,length(long_random_dds_norm$colname),by=10)])

dev.off()

png(file = "heatmap_edge_norm_counts.png", height = 7, width = 10, units = "in", res = 300)
p3 <- ggplot(long_random_edge_norm, aes(x = colname, y = rowname, fill = norm)) +
  geom_tile() +
  labs(x = "Patient ID", y = "Gene ID", tag = "C", fill = "Counts") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    text = element_text(size = 20)
  ) +
  scale_fill_viridis(discrete = FALSE) 
  #scale_x_discrete(breaks=long_random_edge_norm$colname[seq(1,length(long_random_edge_norm$colname),by=10)])
dev.off()


# --------- now with ALDEx2 ------------
# what if we do the same after DESEQ normalization and EDGER
# open the files
aldex_norm <- read.csv(unz("results/bulk/aldex2_clr.csv.zip", "aldex2_clr.csv"))

# do the same plot as before:
# sort the same
aldex_norm <- aldex_norm %>% arrange(X)


# clean data
row.names(aldex_norm) <- aldex_norm$X
aldex_norm <- aldex_norm[, -1]


# Import head of raw bulk counts to latex table
xtable(aldex_norm[c(1:5), c(1:4, 39)])

# heat map normalized counts random sample
# get random cols and rows
random_aldex_norm <- aldex_norm[rows, cols]


# We need a long dataframe to plot the heatmap
long_random_aldex_norm <- random_aldex_norm %>%
  rownames_to_column() %>%
  gather(colname, Counts, -rowname)
head(long_random_aldex_norm)


# Now create the normalized 0-1 column:
long_random_aldex_norm <- long_random_aldex_norm %>% mutate(norm = (Counts - min(Counts)) / (max(Counts) - min(Counts)))

png(file = "heatmap_aldex_norm_counts.png", height = 7, width = 10, units = "in", res = 300)
p4 <- ggplot(long_random_aldex_norm, aes(x = colname, y = rowname, fill = norm)) +
  geom_tile() +
  labs(x = "Patient ID", y = "Gene ID", tag = "D", fill = "Counts") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    text = element_text(size = 20)
  ) +
  scale_fill_viridis(discrete = FALSE) #+
#scale_x_discrete(breaks=long_random_dds_norm$colname[seq(1,length(long_random_dds_norm$colname),by=10)])
dev.off()


require(gridExtra)
png("foo_counts.png", height = 2100, width = 2100)
grid.arrange(p1, p2, p3, p4, ncol = 2, nrow=2)
dev.off()


require("egg")
png("foo_counts_2.png", height = 2100, width = 2100)
ggarrange(p1, p2, p3, p4, ncol = 2, nrow=2)
dev.off()


