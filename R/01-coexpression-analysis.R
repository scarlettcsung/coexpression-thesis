source("R/packages.R")
source("R/coexpression-functions.R")

# ================= DATA SETUP =================

path <- "data/k12_counts.csv" # File Path
raw_counts <- read.csv(path, header=TRUE, row.names = 1) # Load .csv
raw_counts <- raw_counts[order(colnames(raw_counts))]

# ================= CO-EXPRESSION ANALYSIS =================

# Filtered Counts
counts <- filterLowExpression(raw_counts) #Rows: Gene name, Columns: Conditions
counts_filtered <- counts[rowSums(counts) > 0, colSums(counts) > 0]

k12_coexpr <- coexpr_analysis(counts_filtered,2)
matrices <- k12_coexpr$matrix
genepairs <- k12_coexpr$genepairs

plot_all_methods_by_bins(counts_filtered,2)

saveRDS(matrices,file="results/k12_matrices.rds")
saveRDS(genepairs,file="results/k12_genepairs.rds")




