# Install and load package "librarian" to install and load other packages 
# efficiently
if (!requireNamespace("librarian", quietly = TRUE)) 
  install.packages("librarian")
library(librarian)

# Install and load package "COBAIN" (older version of COCOLOCUS)
if (!requireNamespace("COBAIN", quietly = TRUE)) {
  install.packages(
    "dependencies/COBAIN.tar.gz",
    repos = NULL,
    type = "source"
  )
}
library(COBAIN)

# Install and load package "BiocManager" which is required for many other 
# bioinformatics R-packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
library(BiocManager)

# Install and load other required packages
librarian::shelf(
  # Co-expression Analysis
  compositions,
  zCompositions,
  tpq/propr,
  WGCNA,
  # Data Formatting
  dplyr,
  tibble,
  data.table,
  tidyverse,
  # Plotting
  ggplot2,
  ggridges,
  patchwork
)


library(GO.db)
library(AnnotationDbi)
library(ontologyIndex)
library(bench)
