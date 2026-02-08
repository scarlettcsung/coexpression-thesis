source("R/packages.R")
source("R/permutation_testing_functions.R")

df <- readRDS("results/k12_GOmatches.rds")

n_permutations <- 5000
n_genepairs <- floor(nrow(df$Original) * 0.001)
mean_allpairs <- df$Original[nrow(df$Original),"GO_match_cumulative"]/100


background_means <- generate_background_means(allgenepairs_df = df$Original,
                                              n_permutations = n_permutations,
                                              n_genepairs = n_genepairs,
                                              col = "GO_overlap",
                                              seed = 41)

saveRDS(background_means,"evaluation/k12_background_means_5000perms.rds")

p_vals <- sapply(names(df),function(nm) {
  df_method <- df[[nm]]
  mean_subset <- df_method[n_genepairs,"GO_match_cumulative"]
  calc_background_p(background_means = background_means,
                    mean_allpairs = mean_allpairs,
                    mean_subset = mean_subset,
                    n_permutations = n_permutations)
})
  
saveRDS(p_vals,"evaluation/k12_pvals_5000perms.rds")  
