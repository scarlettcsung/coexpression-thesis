source("R/packages.R")
source("R/plotting-functions.R")

# Re-order matrices 
matrices_burk <- readRDS("results_ignore/burkholderia_matrices.rds")
burk_sorted_names <- readRDS("reference/derived/burkholderia_genenames_sorted.rds")
burk_ordered_matrices <- lapply(matrices_burk, function(m) {
  reorder_matrix(m, burk_sorted_names)
})

matrices_k12 <- readRDS("results_ignore/k12_matrices.rds")
k12_sorted_names <- readRDS("reference/derived/k12_genenames_sorted.rds")
k12_ordered_matrices <- lapply(matrices_k12, function(m) {
  reorder_matrix(m, k12_sorted_names)
})

matrices_bac <- readRDS("results_ignore/bacillus_matrices.rds")
bac_sorted_names <- readRDS("reference/derived/bacillus_genenames_sorted.rds")
bac_ordered_matrices <- lapply(matrices_bac, function(m) {
  reorder_matrix(m, bac_sorted_names)
})

# Extract correlations by bins
binned_burk <- extract_bin_correlations(burk_ordered_matrices)
binned_k12 <- extract_bin_correlations(k12_ordered_matrices)
binned_bac <- extract_bin_correlations(bac_ordered_matrices)

method_order <- c("Original", "MeanResiduals", "PC1Residuals", 
              "CLR", "propr","propr_z")
dataset_order <- c("B_pseudomallei","E_coli","B_subtillis")

binned <- list(B_pseudomallei = binned_burk,
               E_coli = binned_k12,
               B_subtillis = binned_bac)

binned_combined <- dplyr::bind_rows(binned, .id = "dataset")
saveRDS(binned_combined, file = "evaluation/eval_results/binned_combined.rds")


png("evaluation/density_coexpression.png", 
    width = 14, 
    height = 6, 
    units = "in", 
    res = 300)

plot_multidataset_collapsed(binned_combined, 
                            method_order = method_order,
                            dataset_order = dataset_order)

dev.off()

median_summary <- generate_median_table(binned_combined)
View(median_summary)


write.csv(median_summary, 
          "evaluation/coexpression_medians.csv", 
          row.names = FALSE)
