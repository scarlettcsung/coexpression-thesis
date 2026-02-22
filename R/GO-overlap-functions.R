source("R/packages.R")
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