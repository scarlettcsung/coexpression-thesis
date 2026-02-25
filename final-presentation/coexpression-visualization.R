library(tidyverse)

k12_genepairs <- readRDS("results_ignore/k12_genepairs.rds")
k12_counts <- read_csv("data/k12_counts.csv")

genes_to_plot = c("b3517","b1493","b1077","b1938")
subset_df <- k12_counts %>% 
  filter(Geneid %in% genes_to_plot)
