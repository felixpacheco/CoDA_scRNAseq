library("tidyverse")
library("DESeq2")

# Read the data
bulk_counts <- read.csv("/data/bulk/counts_bulk.tsv", sep = "\t")
metadata <- read.csv("/data/bulk/metadata_bulk.tsv", sep = "\t")

# Clean metadata
# remove columns with all NA
metadata <- metadata %>% select_if(~ !all(is.na(.)))
# Rename columns, remove "Sample.Characteristic" and "." at the end of column name
names(metadata) <- gsub("Sample.Characteristic.", "", names(metadata))
names(metadata) <- sub(".$", "", names(metadata))
# Drop unused columns
# Age into only number
metadata <- metadata %>%
  rename("Ru" = "run") %>%
  select(-contains("Ontology.Term")) %>%
  separate(age, c("age", "drop")) %>%
  subset(select = -drop)

# Loading data to DESeq2
