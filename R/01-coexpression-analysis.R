# Author: Scarlet (Tsz Ching) Sung
# Description: 
# Main script for co-expression-analysis.
# Runs co-expression analysis functions.
# Saves filtered pairs, similarity matrices, and sorted gene pairs by similarity

# Load packages
source("R/packages.R")
# Load functions
source("R/coexpression-functions.R")

# ================= DATA SETUP =================
# Here, we load the gene expression matrix file/raw-counts
# Files should be .csv files with locus tag and RNA-seq raw counts
# Locus tags (ex. BPS_RS00005) are row names
# Conditions (ex. Control_1) are column names
# Files should not contain other information

# File path
path <- "data/BURK_counts.csv" 
# Load .csv
raw_counts <- read.csv(path, header=TRUE, row.names = 1) 
# Orders samples by alphabetical order
raw_counts <- raw_counts[order(colnames(raw_counts))] 

# ================= CO-EXPRESSION ANALYSIS =================

# Filtering Counts
counts <- filterLowExpression(raw_counts) 
counts_filtered <- counts[rowSums(counts) > 0, colSums(counts) > 0]

# Optional: saves filtered counts
saveRDS(counts_filtered,file = "results/burkholderia_filtered.rds")

# Optional: loads filtered counts
counts_filtered <- readRDS("results/burkholderia_filtered.rds")

# Runs all co-expression analysis methods
coexpr <- coexpr_analysis(counts = counts_filtered, replicates = 3)

# Assigns matrix and genepairs to variables for saving
matrices <- coexpr$matrix
genepairs <- coexpr$genepairs

# Optional: saves matrices and gene pairs. 
# Not included in repository as data sets are too big for GitHub.
# However, these can be generated consistently via. saved filtered counts.
saveRDS(matrices,file="results_ignore/burkholderia_matrices.rds")
saveRDS(genepairs,file="results_ignore/burkholderia_genepairs.rds")

# ================= FORMATTING FOR CORRELATION DISTRIBUTIONS =================
# This is not strictly co-expression analysis, but since the data are already
# being processed, this will be useful for examining correlation/proportionality
# distributions

# Transposing for sorting expression
counts_t <- t(counts_filtered) 
# Sorting by expression for distribution analysis
counts_sorted <- sort_expression(counts_t) 
# Sorted gene names for distribution analysis
sorted_gene_names <- colnames(counts_sorted)

# Saves sorted gene names
saveRDS(sorted_gene_names,
        file="reference/derived/burkholderia_genenames_sorted.rds")
