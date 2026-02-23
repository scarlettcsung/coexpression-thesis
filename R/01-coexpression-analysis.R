# Author: Scarlet (Tsz Ching) Sung
# Description: 
#   Main script for co-expression-analysis.
#   Runs co-expression analysis functions.
#   Saves filtered pairs, similarity matrices, 
#   and gene pairs sorted by similarity
# Usage: Rscript R/01-coexpression-analysis.R {input_filepath} {dataset_name}
#   Example: Rscript 01-coexpression-analysis.R {data/BURK_counts.csv} {burkholderia}
#   NOTE: You may also run this script in an IDE. In that case,
#         remember to change the input file path and data set name!

# Load packages
source("R/packages.R")
# Load functions
source("R/coexpression-functions.R")

# ================= INPUT AND FILENAME SETUP =================
args <- commandArgs(trailingOnly = TRUE)

# If running from IDE, remember to edit and run this block!
if (length(args) >= 2) {
  # Running from command line
  input_file   <- args[1]
  dataset_name <- args[2]
} else {
  # Running from IDE. Remember to change the input_file and dataset_name!
  input_file   <- "data/BURK_counts.csv"
  dataset_name <- "burkholderia"
}

cat("args:", paste(args, collapse = " "), "\n")

# ================= DATA SETUP =================
# Here, we load the gene expression matrix file/raw-counts
# Files should be .csv files with locus tag and RNA-seq raw counts
# Locus tags (ex. BPS_RS00005) are row names
# Conditions (ex. Control_1) are column names
# Files should not contain other information

# File path
path <- input_file
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

# Optional: loads filtered counts. Commented out, but code is usable.
# counts_filtered <- readRDS(file.path("results", sprintf("%s_filtered.rds", dataset_name)))

# Runs all co-expression analysis methods
coexpr <- coexpr_analysis(counts = counts_filtered, replicates = 3)

# Assigns matrix and genepairs to variables for saving
matrices <- coexpr$matrix
genepairs <- coexpr$genepairs

# Optional: saves matrices and gene pairs. 
# Not included in repository as data sets are too big for GitHub.
# However, these can be generated consistently via. saved filtered counts.
saveRDS(matrices,
        file.path("results_ignore", sprintf("%s_matrices.rds",dataset_name)))
        
saveRDS(genepairs,
        file.path("results_ignore", sprintf("%s_genepairs.rds",dataset_name)))

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
        file.path("reference","derived",
                  sprintf("%s_genenames_sorted.rds",dataset_name)))
