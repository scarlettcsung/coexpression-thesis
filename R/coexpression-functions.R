source("R/packages.R")

# ================= HELPER FUNCTIONS =================

# Sort genes by mean expression
sort_expression <- function(df) {
  dataset_sorted <- df[, names(sort(colMeans(df)))]
  return(dataset_sorted)
}

# Compositions CLR
get_compositions <- function(df) {
  # gene_names <- colnames(df)
  run_clr <- compositions::clr(df)
  clr_proportion <- compositions::clrInv(run_clr)
  # browser()
  # colnames(clr_proportion) <- gene_names
  return(clr_proportion)
}

# Pearson correlation of CLR (WGCNA)
get_pcc <- function(df) {
  cor_mat <- WGCNA::cor(df, method = "pearson" #, use = "pairwise.complete.obs"
  )
  message("correlation calculation complete") # Diagnostic line
  return(cor_mat)
}

# Proportionality (propr)
get_corrs <- function(df) {
  run_propr <- propr(df,
                     metric = "rho",
                     ivar = "clr",
                     alpha = NA,
                     p = 100)
  corr_matrix <- getMatrix(run_propr)
  colnames(corr_matrix) <- colnames(df)
  rownames(corr_matrix) <- colnames(df)
  return(corr_matrix)
}

# Extract gene pairs from matrix
get_genepairs <- function(corr_matrix) {
  # Get lower triangle indices (exclude diagonal)
  idx <- which(lower.tri(corr_matrix, diag = FALSE), arr.ind = TRUE)
  
  # Build data frame of gene pairs and correlations
  corr_df <- data.frame(
    gene1 = rownames(corr_matrix)[idx[, 1]],
    gene2 = colnames(corr_matrix)[idx[, 2]],
    correlation = corr_matrix[idx]
  )
  
  # Sort by descending correlation
  corr_df <- corr_df[order(-corr_df$correlation), ]
  
  return(corr_df)
}

# ================= CORRELATION ANALYSIS FOR BENCHMARKING =================

corrs_Original <- function(counts,replicates) {
  tmm <- logTMM(counts,replicates = replicates)
  sorted <- tmm[order(rowMeans(tmm)), ]
  means <- stats::setNames(rowMeans(sorted), rownames(sorted))
  corr <- corAndSpQN(sorted,means)
  corr <- corr[rownames(counts), rownames(counts), drop = FALSE]
  return(corr)
}

corrs_MeanResiduals <- function(counts,replicates) {
  tmm <- logTMM(counts,replicates = replicates)
  sorted <- tmm[order(rowMeans(tmm)), ]
  lemmed <- meanResiduals(sorted)
  mr <- lemmed[rownames(sorted), ]
  means <- stats::setNames(rowMeans(sorted), rownames(sorted))
  corr <- corAndSpQN(mr,means)
  corr <- corr[rownames(counts), rownames(counts), drop = FALSE]
  return(corr)
}

corrs_PC1Residuals <- function(counts,replicates) {
  tmm <- logTMM(counts,replicates = replicates)
  sorted <- tmm[order(rowMeans(tmm)), ]
  less_pc <- suppressMessages(removePCs(sorted, 1))
  means <- stats::setNames(rowMeans(sorted), rownames(sorted))
  pc1 <- less_pc[rownames(sorted), ]
  corr <- corAndSpQN(pc1,means)
  corr <- corr[rownames(counts), rownames(counts), drop = FALSE]
  return(corr)
}

corrs_CLR <- function(counts) {
  counts_t <- t(counts) 
  counts_sorted <- sort_expression(counts_t) 
  z_counts <- cmultRepl(counts_sorted,method = "CZM",output="p-counts") 
  comp <- get_compositions(z_counts)
  comp_pcc <- get_pcc(comp)
  return(comp_pcc)
}

get_propr <- function(counts) {
  counts_t <- t(counts) 
  counts_sorted <- sort_expression(counts_t) 
  pr <- get_corrs(counts_sorted)
  return(pr)
}

get_zpropr <- function(counts) {
  counts_t <- t(counts) 
  counts_sorted <- sort_expression(counts_t)
  z_counts <- cmultRepl(counts_sorted,method = "CZM",output="p-counts")
  pr_z <- get_corrs(z_counts)
  return(pr_z)
}

# ================= WRAPPER =================
coexpr_analysis <- function(counts, replicates) {
  counts_t <- t(counts) # Transposed for formatting
  counts_sorted <- sort_expression(counts_t) # Sorting genes by mean expression
  z_counts <- cmultRepl(counts_sorted,method = "CZM",output="p-counts") # zCompositions zero imputation
  
  comp <- get_compositions(z_counts)
  comp_pcc <- get_pcc(comp)
  pr <- get_corrs(counts_sorted)
  pr_z <- get_corrs(z_counts)
  
  tmm <- logTMM(counts,replicates = replicates)
  sorted <- tmm[order(rowMeans(tmm)), ]
  lemmed <- meanResiduals(sorted)
  less_pc <- suppressMessages(removePCs(sorted, 1))
  results <- list(Original = sorted, 
                  MeanResiduals = lemmed, 
                  PC1Residuals = less_pc)
  results <- lapply(results, function(x) x[rownames(counts), ])
  
  means <- stats::setNames(rowMeans(sorted), rownames(sorted))
  corr_og <- corAndSpQN(results$Original,means)
  corr_mr <- corAndSpQN(results$MeanResiduals, means)
  corr_pc1 <- corAndSpQN(results$PC1Residuals, means)
  
  matrix <- list()
  matrix$Original = corr_og
  matrix$MeanResiduals = corr_mr
  matrix$PC1Residuals = corr_pc1
  matrix$CLR = comp_pcc
  matrix$propr = pr
  matrix$propr_z = pr_z
  
  genepairs <- list()
  genepairs$Original = get_genepairs(corr_og)
  genepairs$MeanResiduals = get_genepairs(corr_mr)
  genepairs$PC1Residuals = get_genepairs(corr_pc1)
  genepairs$CLR = get_genepairs(comp_pcc)
  genepairs$propr = get_genepairs(pr)
  genepairs$propr_z = get_genepairs(pr_z)
  
  output <- list(matrix = matrix, genepairs = genepairs)
  
  return(output)
}

# ================= DENSITY PLOTS =================


# ================= GO-ID MATCHES =================

filter_process_go <- function(genepairs_df, go_reference) {
  
  # Keep only genes that have at least one GO "process" annotation
  process_genes <- unique(go_reference$locus_tag[go_reference$type == "process"])
  
  # Keep only pairs where BOTH genes have process GO terms
  filtered_pairs <- subset(
    genepairs_df,
    gene1 %in% process_genes & gene2 %in% process_genes
  )
  
  before_nrows <- nrow(genepairs_df)
  filtered_nrows <- nrow(filtered_pairs)
  
  p_remaining = 100 * (filtered_nrows/before_nrows)
  
  message(sprintf("Filtered gene pairs: %d → %d (%.1f%% retained)",
                  before_nrows, filtered_nrows, p_remaining))
  return(filtered_pairs)
}

get_process_go_ids <- function(gene_pairs_df, go_df, 
                               gene1_col = "gene1", gene2_col = "gene2",
                               cor_values = NULL) {
  
  # Keep only process GO terms
  go_process <- go_df[go_df$type == "process", ]
  
  # Precompute GO terms per gene as a named list
  go_lookup <- tapply(go_process$GO_ID, go_process$locus_tag, function(x) paste(unique(x), collapse = ";"))
  
  # Vectorized lookup for each gene
  GO_IDs_gene1 <- go_lookup[gene_pairs_df[[gene1_col]]]
  GO_IDs_gene2 <- go_lookup[gene_pairs_df[[gene2_col]]]
  
  # Replace NULL with NA
  GO_IDs_gene1[sapply(GO_IDs_gene1, is.null)] <- NA
  GO_IDs_gene2[sapply(GO_IDs_gene2, is.null)] <- NA
  
  # Build result
  result <- data.frame(
    gene1 = gene_pairs_df[[gene1_col]],
    gene2 = gene_pairs_df[[gene2_col]],
    GO_IDs_gene1 = unlist(GO_IDs_gene1),
    GO_IDs_gene2 = unlist(GO_IDs_gene2),
    correlation = if (!is.null(cor_values)) cor_values else NA,
    stringsAsFactors = FALSE
  )
  
  return(result)
}

compute_go_overlap <- function(df, sep = ";") {
  
  stopifnot(is.data.frame(df))
  stopifnot(c("GO_IDs_gene1", "GO_IDs_gene2") %in% colnames(df))
  
  # compute binary overlap
  df$GO_overlap <- mapply(function(g1, g2) {
    
    if (is.na(g1) || is.na(g2) || g1 == "" || g2 == "") return(0)
    
    go1 <- strsplit(g1, sep, fixed = TRUE)[[1]]
    go2 <- strsplit(g2, sep, fixed = TRUE)[[1]]
    
    as.integer(length(intersect(go1, go2)) > 0)
    
  }, df$GO_IDs_gene1, df$GO_IDs_gene2)
  
  # running proportion of 1s
  df$GO_match_cumulative <- 100 * (cumsum(df$GO_overlap) / seq_len(nrow(df)))
  message("GO overlap computation complete")
  
  df
}

get_GO_and_corr <- function(row_number, db_list, db_names = NULL) {
  # db_list: list of dataframes
  # db_names: optional names for each database; defaults to 1:length(db_list)
  if (is.null(db_names)) db_names <- paste0("db", seq_along(db_list))
  
  # Build the data.frame
  result <- data.frame(
    GO_match_cumulative = sapply(db_list, function(df) {
      if (row_number > nrow(df)) return(NA)
      df$GO_match_cumulative[row_number]
    }),
    correlation = sapply(db_list, function(df) {
      if (row_number > nrow(df)) return(NA)
      df$correlation[row_number]
    }),
    stringsAsFactors = FALSE
  )
  
  return(result)
}