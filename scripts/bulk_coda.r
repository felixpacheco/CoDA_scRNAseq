# -------------------- Load the libraries------------------------------------
library("tidyverse")
#library("xtable")

# -------------------- Read the data ----------------------------------------
raw_bulk_counts <- read.csv("/home/people/laucom/CoDA_scRNAseq/data/bulk/counts_bulk_patient.tsv", sep = "\t")
metadata_raw <- read.csv("/home/people/laucom/CoDA_scRNAseq/data/bulk/metadata_bulk.tsv", sep = "\t")
# ---------------- Wrangle bulk counts --------------------------------------
bulk_counts <- raw_bulk_counts
bulk_counts[63152, "Ensembl_gene_id"] <- "TC_"
bulk_counts <- bulk_counts[-c(63152), ]
bulk_counts_clean <- bulk_counts[apply(bulk_counts[, -1], 1, function(x) !all(x == 0)), ]

# Import head of raw bulk counts to latex table
#xtable(bulk_counts_clean[c(1:5), c(1:4, 39)])

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

# Wrangle a bit more
row.names(bulk_counts_clean) <- bulk_counts_clean$Ensembl_gene_id
bulk_counts_clean <- bulk_counts_clean[, -1]

# ----------------- ALDEx2 package -------------------------------------------
# ----------------------------------------------------------------------------
# the BioCManager packages need to be loaded always after the data wrangling because they interfere
library("ALDEx2")
conds <- metadata_counts$sample_type
#x.all <- aldex(bulk_counts_clean, conds, mc.samples=128, test="t", effect=TRUE,
               #include.sample.summary=FALSE, denom="all", verbose=FALSE)


# Write to file
#write.csv(x.all, "aldex2_DE.csv", row.names = TRUE)

# clr
x <- aldex.clr(bulk_counts_clean, conds, mc.samples=128, denom="all", verbose=F)

#x_df <- data.frame(x@analysisData)

# now reduce the columns to one
#library(stringr)
#nm1 <- str_replace(names(x_df),"\\..*", "")
#x.new <- data.frame
#x_df[paste0(unique(nm1), "_mean")] <- sapply(split.default(x_df, nm1), rowMeans)
#x.new <- x_df[grep("_mean", names(x_df))] 
#write.csv(x.new, "aldex2_clr.csv", row.names = TRUE)

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

kw <- aldex.kw(x)
kw.effect <- aldex.effect(x, include.sample.summary=TRUE, glm.conds=TRUE, CI=TRUE)
print("kw and kw.effect done")
write.csv(kw, gzfile("aldex2_kw.csv.gz"), row.names = TRUE)
write.csv(kw.effect, gzfile("aldex2_kw_clr.csv.gz"), row.names = TRUE)
print("files saved")


