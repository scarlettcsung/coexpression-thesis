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
  results <- list(Original = sorted, MeanResiduals = lemmed, PC1Residuals = less_pc)
  results <- lapply(results, function(x) x[rownames(counts), ])
  
  means <- stats::setNames(rowMeans(sorted), rownames(sorted))
  corr_og <- corAndSpQN(results$Original,means)
  corr_mr <- corAndSpQN(results$MeanResiduals, means)
  corr_pc1 <- corAndSpQN(results$PC1Residuals, means)
  
  matrix <- list(Original = corr_og, 
                        MeanResiduals = corr_mr, 
                        PC1Residuals = corr_pc1,
                        CLR = comp_pcc,
                        propr = pr,
                        propr_z = pr_z)
  genepairs <- list(Original = get_genepairs(corr_og), 
                  MeanResiduals = get_genepairs(corr_mr), 
                  PC1Residuals = get_genepairs(corr_pc1),
                  CLR = get_genepairs(comp_pcc),
                  propr = get_genepairs(pr),
                  propr_z = get_genepairs(pr_z))
  output <- list(matrix = matrix, genepairs = genepairs)
  
  return(output)
}

# ================= DENSITY PLOTS =================
plot_collapsed_by_subset <- function(all_cors) {
  
  # ---- reshape to long ----
  cors_df <- do.call(rbind, lapply(names(all_cors), function(method) {
    do.call(rbind, lapply(names(all_cors[[method]]), function(subset) {
      data.frame(
        value  = all_cors[[method]][[subset]],
        Method = method,
        Subset = subset,
        stringsAsFactors = FALSE
      )
    }))
  }))
  
  cors_df$Subset <- factor(
    cors_df$Subset,
    levels = c(
      "Bottom 10% to self",
      "Middle 10% to self",
      "Top 10% to self",
      "Bottom 10% to Top 10%"
    )
  )
  
  # ---- plot ----
  ggplot(cors_df, aes(x = value, y = Method, fill = Method)) +
    geom_segment(aes(x = -1, xend = 1, yend = Method), alpha = 0.4) +
    ggridges::geom_density_ridges(
      scale = 0.9,
      quantile_lines = TRUE,
      show.legend = FALSE
    ) +
    facet_wrap(~ Subset, nrow = 1) +
    scale_x_continuous(limits = c(-1, 1)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    coord_flip() +
    theme_minimal(base_size = 7) +
    theme(
      axis.title = element_blank(),
      strip.text = element_text(size = 7),
      axis.text.y = element_text(size = 6)
    )
}

plot_all_methods_by_bins <- function(counts, replicates) {
  
  # ---------------- Helper functions ----------------
  ..sort_bins <- function(expr_matrix) {
    dec <- floor(nrow(expr_matrix)/10)
    bins <- list(
      low  = expr_matrix[seq(1, length.out = dec), ],
      mid  = expr_matrix[seq(from = (nrow(expr_matrix) - dec)/2, length.out = dec), ],
      top  = expr_matrix[seq(to = nrow(expr_matrix), length.out = dec), ]
    )
    combined <- rbind(bins$low, bins$top)
    list(bins = bins, cross = combined)
  }
  
  ..gatherCors <- function(expr_matrix) {
    bin_data <- ..sort_bins(expr_matrix)
    bins <- bin_data$bins
    cross <- bin_data$cross
    
    stats::setNames(mapply(function(x, y) {
      mat <- WGCNA::cor(t(bins[[x]]), t(bins[[y]]))
      mat[upper.tri(mat)]
    }, x = c("low","mid","top","low"), y = c("low","mid","top","top"), SIMPLIFY = F),
    c("Bottom 10% to self", "Middle 10% to self", "Top 10% to self", "Bottom 10% to Top 10%"))
  }
  
  ..gatherFullCors <- function(corr_matrix) {
    n <- nrow(corr_matrix)
    dec <- floor(n / 10)
    indices <- list(
      low = 1:dec,
      mid = floor((n - dec)/2 + 1):(floor((n - dec)/2 + 1) + dec - 1),
      top = (n - dec + 1):n
    )
    list(
      "Bottom 10% to self" = corr_matrix[indices$low, indices$low] |> as.vector(),
      "Middle 10% to self" = corr_matrix[indices$mid, indices$mid] |> as.vector(),
      "Top 10% to self" = corr_matrix[indices$top, indices$top] |> as.vector(),
      "Bottom 10% to Top 10%" = corr_matrix[indices$low, indices$top] |> as.vector()
    )
  }
  
  ..denPlot <- function(cor_list) {
    ..corToRGB <- function(cor) {
      cor_colors <- c("#67001F","#B2182B","#D6604D","#F4A582","#FDDBC7","#FFFFFF",
                      "#D1E5F0","#92C5DE","#4393C3","#2166AC","#053061")
      to_color <- grDevices::colorRamp(cor_colors)
      grDevices::rgb(to_color((cor+1)/2), maxColorValue = 255)
    }
    medians <- lapply(cor_list, median)
    fills <- unlist(lapply(medians, ..corToRGB))
    
    ggplot(utils::stack(cor_list), aes(x = values, y = ind, fill = ind)) +
      geom_segment(mapping = aes(x=-1, xend=1, y=ind, yend=ind)) +
      ggridges::geom_density_ridges(scale=0.95, quantile_lines=TRUE, show.legend = FALSE) +
      scale_fill_manual(values=fills) +
      ggridges::theme_ridges(center_axis_labels = TRUE) +
      theme(axis.text.x=element_text(size = 1, hjust=0),
            axis.text.y = element_text(size = 1),
            axis.title.x=element_blank(),
            axis.title.y=element_blank()) +
      scale_x_continuous(limits=c(-1,1)) +
      scale_y_discrete(expand = expansion(mult=c(0.03,0.05))) +
      geom_vline(xintercept = 0, linetype="dashed") +
      coord_flip()
  }
  
  # ---------------- TMM-based ----------------
  tmm <- logTMM(counts, replicates = replicates)
  sorted <- tmm[order(rowMeans(tmm)), ]
  lemmed <- meanResiduals(sorted)
  less_pc <- suppressMessages(removePCs(sorted, 1))
  
  tmm_results <- list(
    Original = sorted,
    MeanResiduals = lemmed,
    PC1Residuals = less_pc
  )
  
  tmm_cors <- lapply(tmm_results, ..gatherCors)
  
  # ---------------- CLR / propr / propr_z ----------------
  counts_t <- t(counts)
  counts_sorted <- sort_expression(counts_t)
  z_counts <- cmultRepl(counts_sorted, method="CZM", output="p-counts")
  
  CLR <- get_pcc(get_compositions(z_counts))
  propr_mat <- get_corrs(counts_sorted)
  propr_z_mat <- get_corrs(z_counts)
  
  extra_cors <- list(
    CLR = ..gatherFullCors(CLR),
    propr = ..gatherFullCors(propr_mat),
    propr_z = ..gatherFullCors(propr_z_mat)
  )
  
  # ---------------- Combine all ----------------
  all_cors <- c(tmm_cors, extra_cors)
  
  # ---------------- Plot everything ----------------
  print(plot_collapsed_by_subset(all_cors))
}

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