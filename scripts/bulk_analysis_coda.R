# -------------------- Load the libraries------------------------------------
library("tidyverse")

# -------------------- Read the data ----------------------------------------
raw_bulk_counts <- read.csv("data/bulk/counts_bulk_patient.tsv", sep = "\t")
metadata <- read.csv("data/bulk/metadata_bulk.tsv", sep = "\t")
# ---------------- Wrangle bulk counts --------------------------------------
bulk_counts <- raw_bulk_counts %>%
  column_to_rownames(var = "Ensembl_gene_id")
# -------------------- Clean metadata ------------------q---------------------
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

# ----------------- ALDEx2 package -------------------------------------------
# ----------------------------------------------------------------------------
# the BioCManager packages need to be loaded always after the data wrangling because they interfere
library("ALDEx2")
conds <- metadata_counts$sample_type
x.all <- aldex(bulk_counts, conds, mc.samples=128, test="t", effect=TRUE,
               include.sample.summary=FALSE, denom="all", verbose=FALSE)

# MC 128 because the t-test will be accurate


