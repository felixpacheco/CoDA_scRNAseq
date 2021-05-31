library("tidyverse")
library("grDevices")

# ----------------- READ THE METADATA --------------------
metadata <- read.csv("/Users/laurasansc/Desktop/BioinformaticsSystemsBiology/CoDA_special_course/data/FACS/annotations_FACS.csv")

# clean metadata
types_count <- metadata %>%
  group_by(tissue) %>%
  count(cell_ontology_class)

# make a plot
type.graph <- function(df2, na.rm = TRUE, ...) {
  
  # Create list of locations
  tissue_list <-unique(df2$tissue)
  
  # Create a for loop to produce ggpot plots
  for (i in seq_along(tissue_list)) {
    filename <- paste(tissue_list[i], "png", sep = ".")
    png(file = filename, height = 7, width = 10, units = "in", res = 150)
    # create a plot for each loc in df
    plot<-
      ggplot(subset(df2, df2$tissue == tissue_list[i]),
             aes(x = cell_ontology_class, 
                 y = n)) +
      geom_bar(stat = "identity") +
      ggtitle(tissue_list[i]) #+
      facet_wrap(i)
    #windows()
    print(plot)
    dev.off()
  }
}

type.graph(types_count)


# We decided to go for heart - fibroblast / fat - mesenchymal stem cell of adipose
# We need to clean the metadata and merge the counts 
counts_heart <- read.csv("/Users/laurasansc/Desktop/BioinformaticsSystemsBiology/CoDA_special_course/data/FACS/Heart-counts.csv")
counts_fat <- read.csv("/Users/laurasansc/Desktop/BioinformaticsSystemsBiology/CoDA_special_course/data/FACS/Fat-counts.csv")


# clean metadata, only keep cells of interest
metadata <- metadata[metadata$tissue == "Heart" | metadata$tissue == "Fat", ] 
metadata <- metadata[metadata$cell_ontology_class == "fibroblast" | metadata$cell_ontology_class == "mesenchymal stem cell of adipose", ] 

# get list of cells of interest
cell_ids <- metadata$cell

# make genes the index
row.names(counts_fat) <- counts_fat$gene
counts_fat <- counts_fat[, -1]

row.names(counts_heart) <- counts_heart$gene
counts_heart <- counts_heart[, -1]

# remove columns that are not our cells of interest
counts_heart <- counts_heart[, colnames(counts_heart) %in% cell_ids]
counts_fat <- counts_fat[, colnames(counts_fat) %in% cell_ids]

# merge
counts <- merge(counts_heart, counts_fat, by=0, all= TRUE)

# save CSV files
write.csv(metadata, gzfile("metadata_sc.csv.gz"), row.names = TRUE)
write.csv(counts, gzfile("counts_sc.csv.gz"), row.names = TRUE)












