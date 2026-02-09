source("R/packages.R")
source("R/coexpression-functions.R")

matrices_burk <- readRDS("results/burkholderia_matrices.rds")
burk_sorted_names <- readRDS("reference/derived/burkholderia_genenames_sorted.rds")
burk_ordered_matrices <- lapply(matrices_burk, function(m) {
  reorder_matrix(m, burk_sorted_names)
})

matrices_k12 <- readRDS("results/k12_matrices.rds")
k12_sorted_names <- readRDS("reference/derived/k12_genenames_sorted.rds")
k12_ordered_matrices <- lapply(matrices_k12, function(m) {
  reorder_matrix(m, k12_sorted_names)
})

matrices_bac <- readRDS("results/bacillus_matrices.rds")
bac_sorted_names <- readRDS("reference/derived/bacillus_genenames_sorted.rds")
bac_ordered_matrices <- lapply(matrices_bac, function(m) {
  reorder_matrix(m, bac_sorted_names)
})

binned_burk <- extract_bin_correlations(burk_ordered_matrices)
binned_k12 <- extract_bin_correlations(k12_ordered_matrices)
binned_bac <- extract_bin_correlations(bac_ordered_matrices)

my_order <- c("Original", "MeanResiduals", "PC1Residuals", "CLR", "propr")

binned <- list(Bpseudomallei = binned_burk,
               Ecoli = binned_k12,
               Bsubtillis = binned_bac)

binned_combined <- dplyr::bind_rows(binned, .id = "dataset")

plot_multidataset_collapsed(binned_combined, method_order = my_order)
