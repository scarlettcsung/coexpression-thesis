source("R/packages.R")
source("R/KOID-overlap-functions.R")

# ================= LOAD DATA =================

genepairs <- readRDS("results_ignore/k12_genepairs.rds")
ko_reference <- read_csv("reference/derived/k12_operon_reference.csv")

# ================= KOID MATCHES =================

# Map operons and notates shared operons
shared_ko <- lapply(genepairs, compute_shared_ko,
                             gene_to_ko = ko_reference)

# Save file for faster use in the future
saveRDS(shared_ko,file="results_ignore/k12_shared_kos.rds")

# Removing pairs where one or both genes are not mapped to operon
filtered_shared_ko <- lapply(shared_ko, filter_missing_ko)

# Adding the cumulative percentages column
ko_matches <- lapply(filtered_shared_ko, add_cumulative_pct)


saveRDS(ko_matches,file="results_ignore/k12_komatches.rds")
