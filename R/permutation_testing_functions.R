source("R/packages.R")

# ================= HELPER FUNCTIONS =================

# ================= PERMUTATION TESTING =================

generate_background_means <- function(allgenepairs_df,
                         n_permutations,
                         n_genepairs,
                         col = "GO_overlap",
                         seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  means <- numeric(n_permutations)
  
  values <- allgenepairs_df[[col]]
  
  for (i in seq_len(n_permutations)) {
    sample_vals <- sample(values, size = n_genepairs, replace = FALSE)
    means[i] <- mean(sample_vals)
  }
  
  return(means)

}

calc_background_p <- function(background_means,
                              mean_allpairs,
                              mean_subset,
                              n_permutations = NULL) {
  
  p_delta <- mean_subset - mean_allpairs
  p_value <- mean((background_means - mean_allpairs) >= p_delta)
  
  if (is.null(n_permutations)) {
    n_permutations <- length(background_means)
  }
  
  if (p_value == 0) {
    message(sprintf(
      "p-value is < %.4g (1/%d permutations)",
      1/n_permutations,
      n_permutations
    ))
  }
  
  return(p_value)
}
