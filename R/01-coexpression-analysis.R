source("R/packages.R")
source("R/coexpression-functions.R")

# ================= DATA SETUP =================

path <- "data/BURK_counts.csv" # File Path
raw_counts <- read.csv(path, header=TRUE, row.names = 1) # Load .csv
raw_counts <- raw_counts[order(colnames(raw_counts))]

# ================= CO-EXPRESSION ANALYSIS =================

# Filtered Counts
counts <- filterLowExpression(raw_counts) #Rows: Gene name, Columns: Conditions
counts_filtered <- counts[rowSums(counts) > 0, colSums(counts) > 0]

saveRDS(counts_filtered,file = "results/bacillus_filtered.rds")

counts_t <- t(counts_filtered)
counts_sorted <- sort_expression(counts_t)
sorted_gene_names <- colnames(counts_sorted)
saveRDS(sorted_gene_names,file="reference/derived/bacillus_genenames_sorted.rds")

coexpr <- coexpr_analysis(counts_filtered,3)
matrices <- coexpr$matrix
genepairs <- coexpr$genepairs

saveRDS(matrices,file="results/bacillus_matrices.rds")
saveRDS(genepairs,file="results/bacillus_genepairs.rds")
