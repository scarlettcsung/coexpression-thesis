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

saveRDS(counts_filtered,file = "results/burkholderia_filtered.rds")

burk_coexpr <- coexpr_analysis(counts_filtered,3)
matrices <- burk_coexpr$matrix
genepairs <- burk_coexpr$genepairs

saveRDS(matrices,file="results/burkholderia_matrices.rds")
saveRDS(genepairs,file="results/burkholderia_genepairs.rds")

plot_all_methods_by_bins(counts_filtered,3)



