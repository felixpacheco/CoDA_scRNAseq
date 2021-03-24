# -------------------- Load the libraries------------------------------------
library("tidyverse")

# -------------------- Read the data ----------------------------------------
raw_bulk_counts <- read.csv("data/bulk/counts_bulk.tsv", sep = "\t")
metadata <- read.csv("data/bulk/metadata_bulk.tsv", sep = "\t")
# ---------------- Wrangle bulk counts --------------------------------------
bulk_counts <- raw_bulk_counts %>%
  subset(select = -Gene.Name) %>%
  column_to_rownames(var = "Gene.ID")
# -------------------- Clean metadata ------------------q---------------------
# remove columns with all NA
metadata <- metadata %>% select_if(~ !all(is.na(.)))
# Rename columns, remove "Sample.Characteristic" and "." at the end of column name
names(metadata) <- gsub("Sample.Characteristic.", "", names(metadata))
names(metadata) <- sub(".$", "", names(metadata))
# Drop unused columns
# Age into only number
metadata <- metadata %>%
  rename(c("run" = "Ru", "sample_type" = "sampling.site")) %>%
  select(-contains("Ontology.Term")) %>%
  separate(age, c("age", "drop")) %>%
  subset(select = -drop)
# ----------------------------------------------------------------------------
# ----------------- ALDEx2 package -------------------------------------------
# ----------------------------------------------------------------------------
library("ALDEx2")
conds <- metadata$sample_type
x.all <- aldex(bulk_counts, conds, mc.samples=16, test="t", effect=TRUE,
               include.sample.summary=FALSE, denom="all", verbose=FALSE)

