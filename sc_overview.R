library("tidyverse")
library("grDevices")
metadata <- read.csv("/Users/laurasansc/Desktop/BioinformaticsSystemsBiology/CoDA_special_course/data/FACS/annotations_FACS.csv")
counts_facs <- read.csv("/Users/laurasansc/Desktop/BioinformaticsSystemsBiology/CoDA_special_course/data/FACS/Brain_Microglia-counts.csv")

# clean metadata
types_count <- metadata %>%
  group_by(tissue) %>%
  count(cell_ontology_class)

# make a plot
tissues <- unique(types_count$tissue)
for (t in tissues)
{
  filename <- paste(t, "png", sep = ".")
  png(file = filename, height = 7, width = 10, units = "in", res = 150)
  types_count %>%
    group_by(tissue) %>%
    nest() %>%
    mutate(plot = map(
      data,
      ~ ggplot(., aes(x = cell_ontology_class, y = n)) +
        geom_bar(stat = "identity") +
        ggtitle(.)
    ))
  dev.off()
}

ggplot(types_count, aes(x = cell_ontology_class, y = n)) +
  geom_bar(stat = "identity") +
  facet_wrap(vars(tissue))


by_type <- types_count %>%
  group_by(tissue) %>%
  nest() %>%
  mutate(plot = map(data,~ ggplot(., aes(x = cell_ontology_class, y = n)) +
      geom_bar(stat = "identity") +
      ggtitle(.) +
      facet_wrap(vars(tissue))
  ))



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





















