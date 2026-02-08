source("R/packages.R")
source("R/coexpression-functions.R")

library(bench)

path <- "data/bacillus_counts.csv" # File Path
raw_counts <- read.csv(path, header=TRUE, row.names = 1) # Load .csv
raw_counts <- raw_counts[order(colnames(raw_counts))]

# ================= CO-EXPRESSION ANALYSIS =================

# Filtered Counts
counts <- filterLowExpression(raw_counts) #Rows: Gene name, Columns: Conditions
counts_filtered <- counts[rowSums(counts) > 0, colSums(counts) > 0]
replicates = 3

# Benchmark the top 3 correlation methods
bm <- bench::mark(
  Original       = corrs_Original(counts, replicates),
  MeanResiduals  = corrs_MeanResiduals(counts, replicates),
  PC1Residuals   = corrs_PC1Residuals(counts, replicates),
  CLR            = corrs_CLR(counts),
  propr          = get_propr(counts),
  propr_z        = get_zpropr(counts),
  iterations = 10,       # repeat each function 10 times
  check = FALSE,          # skip automatic equality checks
  gc = TRUE
)

# View results
print(bm)

bm %>%
  dplyr::select(expression, time) %>%
  group_by(expression) %>%
  summarise(
    mean_time = mean(as.numeric(as_bench_time(unlist(time)), units = "ms")),
    sd_time   = sd(as.numeric(as_bench_time(unlist(time)), units = "ms"))
  )
