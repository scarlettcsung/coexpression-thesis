source("R/packages.R")

# ================= PERMUTATION TESTING =================
#' @description
#' Generates background match means from randomized subsets
#' 
#' @param allgenepairs_df A data frame with gene pairs and match notations (0 or 1)
#' @param n_permutations Numeric, number of permutations
#' @param col String, column name
#' @param seed numeric integer for seed. Default NULL.
#'
#' @return A numeric vector of background means, length of n_permutations
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

#' @description
#' Calculates p-value with permutation testing against background means
#' 
#' @param background_means A numeric vector of background means
#' @param mean_allpairs Numeric, mean overlap of all gene pairs
#' @param mean_subset Numeric, mean overlap of subset of interest
#' @param n_permutations Numeric, number of permutations. Default NULL where 
#'    permutations will be automatically calculated
#'
#' @return Numeric of p-value. If p-value calculated is 0, a message will be
#'    printed with estimated p-value. 
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
