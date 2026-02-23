# Load packages
source("R/packages.R")

# ================= GO-ID MATCHES =================
# GPT 5.1 was used throughout this section for generating and debugging
# functions, which were reviewed and edited manually. 
# Logic was given to GPT for implementation in R.

#' @description
#' Removes gene pairs where either or both gene(s) do not annotate to a 
#' "biological process" GO-term.
#' 
#' @param genepairs_df A data frame describing similarity/"correlation" values
#'    of gene pairs. It has three required columns, described below.
#'    - gene 1: gene name or locus tag of gene 1 of gene pair
#'    - gene 2: gene name or locus tag of gene 2 of gene pair
#'    - correlation: similarity measure of gene pair
#'    
#' @param go_reference A data frame with GO-term information of annotated genes.
#'    Required columns are described below.
#'    - locus tag: locus tag or gene name of gene
#'    - type: type of GO term, e.g. "process" for Biological Process, "function"
#'            for Molecular Function, "component" for Cellular Component
#'    - GO_ID: GO-ID in the format of GO:1234567 (ex.	GO:0008652)
#'    - GO_term: GO-term
#'
#' @return A data frame where only gene pairs with both genes annotated to at 
#'    least one GO-term are retained. Columns same as genepairs_df input.
#' @return Message with statistics regarding gene pairs filtering.
filter_process_go <- function(genepairs_df, go_reference) {
  
  # Keep only genes that have at least one GO "process" annotation
  process_genes <- unique(go_reference$locus_tag[go_reference$type == "process"])
  
  # Keep only pairs where BOTH genes have process GO terms
  filtered_pairs <- subset(
    genepairs_df,
    gene1 %in% process_genes & gene2 %in% process_genes
  )
  
  # Calculating percentage remaining
  before_nrows <- nrow(genepairs_df)
  filtered_nrows <- nrow(filtered_pairs)
  p_remaining = 100 * (filtered_nrows/before_nrows)
  
  message(sprintf("Filtered gene pairs: %d → %d (%.1f%% retained)",
                  before_nrows, filtered_nrows, p_remaining))
  return(filtered_pairs)
}

#' @description
#' Adds "Biological Process" GO-IDs to gene pairs. 
#' 
#' @param genepairs_df A data frame describing similarity/"correlation" values
#'    of gene pairs. It has three required columns, described below.
#'    - gene 1: gene name or locus tag of gene 1 of gene pair
#'    - gene 2: gene name or locus tag of gene 2 of gene pair
#'    - correlation: similarity measure of gene pair
#' @param go_df A data frame with GO-term information of annotated genes.
#'    Required columns are described below.
#'    - locus tag: locus tag or gene name of gene
#'    - type: type of GO term, e.g. "process" for Biological Process, "function"
#'            for Molecular Function, "component" for Cellular Component
#'    - GO_ID: GO-ID in the format of GO:1234567 (ex.	GO:0008652)
#'    - GO_term: GO-term
#' @param gene1_col A string, name of gene 1 column. Default "gene1"
#' @param gene2_col A string, name of gene 2 column. Default "gene2"
#' @param cor_values A string, name of correlation/similarity measure value.
#'    Default NULL - no correlation column will be included. 
#'
#' @return A data frame where genes of gene pairs are annotated to its GO-IDs. 
#'    Columns are described below:
#'    - gene 1: gene name or locus tag of gene 1 of gene pair
#'    - gene 2: gene name or locus tag of gene 2 of gene pair
#'    - GO_IDs_gene1: GO-IDs annotated to gene 1
#'    - GO_IDs_gene2: GO-IDs annotated to gene 2
#'    - correlation: similarity measure of gene pair
get_process_go_ids <- function(gene_pairs_df, go_df, 
                               gene1_col = "gene1", gene2_col = "gene2",
                               cor_values = NULL) {
  
  # Keep only process GO terms. This is redundant to filtered pairs, but could
  # still be useful if genes were not pre-filtered.
  go_process <- go_df[go_df$type == "process", ]
  
  # Precompute GO terms per gene as a named list, to speed up annotations
  go_lookup <- tapply(go_process$GO_ID, go_process$locus_tag, 
                      function(x) paste(unique(x), collapse = ";"))
  
  # Vectorized lookup for each gene
  GO_IDs_gene1 <- go_lookup[gene_pairs_df[[gene1_col]]]
  GO_IDs_gene2 <- go_lookup[gene_pairs_df[[gene2_col]]]
  
  # Replace NULL with NA
  GO_IDs_gene1[sapply(GO_IDs_gene1, is.null)] <- NA
  GO_IDs_gene2[sapply(GO_IDs_gene2, is.null)] <- NA
  
  # Build result data frame
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

#' @description
#' Adds binary notation to notate GO overlap and calculates cumulative overlaps
#' 
#' @param df A data frame where genes of gene pairs are matched to its annotated
#'    GO-IDs. Columns are described below:
#'    - gene1: gene name or locus tag of gene 1 of gene pair
#'    - gene2: gene name or locus tag of gene 2 of gene pair
#'    - GO_IDs_gene1: GO-IDs annotated to gene 1
#'    - GO_IDs_gene2: GO-IDs annotated to gene 2
#'    - correlation: similarity measure of gene pair
#' @param sep A string, separator of GO-IDs in GO_IDs_gene1 and GO_IDs_gene2 of
#'    df. eg. ";" or ",". Default ";".
#'
#' @return A data frame that is identical to "df" input but includes two
#'    additional columns.
#'    - GO_overlap: binary notation 0 for no overlap, 1 for overlap
#'    - GO_match_cumulative: cumulative overlap up to row in percentage
compute_go_overlap <- function(df, sep = ";") {
  
  # Checks for data frame and columns
  stopifnot(is.data.frame(df))
  stopifnot(c("GO_IDs_gene1", "GO_IDs_gene2") %in% colnames(df))
  
  # compute binary overlap
  df$GO_overlap <- mapply(function(g1, g2) {
    
    # check for missing values
    if (is.na(g1) || is.na(g2) || g1 == "" || g2 == "") return(0)
    
    # Splits GO-IDs
    go1 <- strsplit(g1, sep, fixed = TRUE)[[1]]
    go2 <- strsplit(g2, sep, fixed = TRUE)[[1]]
    
    # as.integer:TRUE --> 1, FALSE --> 0
    # so if there are intersects...
    as.integer(length(intersect(go1, go2)) > 0) 
    
  }, df$GO_IDs_gene1, df$GO_IDs_gene2) # working on these columns
  
  # Cumulate % of 1s in GO_overlap up to row number
  df$GO_match_cumulative <- 100 * (cumsum(df$GO_overlap) / seq_len(nrow(df)))
  
  message("GO overlap computation complete")
  return(df)
}

