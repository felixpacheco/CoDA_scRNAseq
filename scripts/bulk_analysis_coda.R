# -------------------- Load the libraries------------------------------------
library("tidyverse")

# -------------------- Read the data ----------------------------------------
raw_bulk_counts <- read.csv("data/bulk/counts_bulk_patient.tsv", sep = "\t")
xtable(raw_bulk_counts[c(1:5), c(1:4, 39)])
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

xtable(metadata_counts[c(1:5),])
# Wrangle a bit more
row.names(bulk_counts_clean) <- bulk_counts_clean$Ensembl_gene_id
bulk_counts_clean <- bulk_counts_clean[, -1]

# ----------------- ALDEx2 package -------------------------------------------
# ----------------------------------------------------------------------------
# the BioCManager packages need to be loaded always after the data wrangling because they interfere
library("ALDEx2")
conds <- metadata_counts$sample_type
x.all <- aldex(bulk_counts, conds, mc.samples=128, test="t", effect=TRUE,
               include.sample.summary=FALSE, denom="all", verbose=FALSE)