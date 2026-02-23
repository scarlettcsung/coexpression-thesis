# Load packages
source("R/packages.R")

# ================= KOID MATCHES =================
# GPT 5.2 was used throughout this section for speeding up generating and 
# debugging functions, all of which were reviewed and edited manually. 
# Logic was given to GPT for implementation in R.
# No unit testing, but outputs were manually inspected against operon database
# NOTE: it might be better to make these functions more consistent with GO-ID
#       overlaps. to fix at a later time. 

#' @description
#' Annotates operon IDs to genes, adds binary notation to notate GO overlap and 
#' adds shared operons between gene pairs
#' @note
#' to-do: make this function logic more consistent with GO-ID matches
#' 
#' @param df A data frame describing similarity/"correlation" values
#'    of gene pairs. It has three required columns, described below.
#'    - gene1: gene name or locus tag of gene 1 of gene pair
#'    - gene2: gene name or locus tag of gene 2 of gene pair
#'    - correlation: similarity measure of gene pair
#' @param gene_to_ko A data frame with operon information of annotated genes.
#'    Required columns are described below.
#'    - gene: locus tag or gene name of gene
#'    - ko: KOID in the format of KO12345 (ex.	KO00009)
#' @param gene1_col A string, name of gene 1 column. Default "gene1"
#' @param gene2_col A string, name of gene 2 column. Default "gene2"
#' @param sep A string, separator to be used of KOIDs in ko_gene1 and ko_gene2 
#'    of df. eg. ";" or ",". Default ",".
#'
#' @return A data frame where gene pairs are annotated to operons. Gene pairs 
#'    with shared operons are denoted, alongside a list of their shared operons. 
#'    Columns are described below:
#'    - gene1: gene name or locus tag of gene 1 of gene pair
#'    - gene2: gene name or locus tag of gene 2 of gene pair
#'    - correlation: similarity measure of gene pair
#'    - ko_gene1: operon ID or list of operon IDs annotated to gene 1. 
#'                NULL if none.
#'    - ko_gene2: operon ID or list of operon IDs annotated to gene 2
#'                NULL if none.
#'    - shared ko: operon ID or list of operon IDs annotated to both genes.
#'                character(0) if none.
#'    - shared_operon: binary notation 0 for no overlap, 1 for overlap
compute_shared_ko <- function(df, gene_to_ko,
                              gene1_col = "gene1",
                              gene2_col = "gene2",
                              sep = ",") {
  
  # build map: gene -> character vector of KOs
  ko_map <- setNames(
    strsplit(gene_to_ko$ko, sep),
    gene_to_ko$gene
  )
  
  # KO lists for each gene
  ko1 <- ko_map[df[[gene1_col]]]
  ko2 <- ko_map[df[[gene2_col]]]
  
  # shared KO (operon) per row
  shared_ko <- mapply(
    function(x, y) {
      if (is.null(x) || is.null(y)) return(character(0))
      intersect(x, y)
    },
    ko1, ko2,
    SIMPLIFY = FALSE
  )
  
  # binary indicator
  shared_binary <- as.integer(lengths(shared_ko) > 0)
  
  # append columns
  df$ko_gene1 <- ko1
  df$ko_gene2 <- ko2
  df$shared_ko <- shared_ko
  df$shared_operon <- shared_binary
  
  message("Shared KO computation complete")
  df
}

#' @description
#' Filters gene pair dataframe made by function compute_shared_ko by removing 
#' gene pairs where one gene or both genes are not annotated to an operon.
#' 
#' @param df A data frame where gene pairs are annotated to operons. Gene pairs 
#'    with shared operons are denoted, alongside a list of their shared operons. 
#'    Columns are described below:
#'    - gene1: gene name or locus tag of gene 1 of gene pair
#'    - gene2: gene name or locus tag of gene 2 of gene pair
#'    - correlation: similarity measure of gene pair
#'    - ko_gene1: operon ID or list of operon IDs annotated to gene 1. 
#'                NULL if none.
#'    - ko_gene2: operon ID or list of operon IDs annotated to gene 2
#'                NULL if none.
#'    - shared ko: operon ID or list of operon IDs annotated to both genes.
#'                character(0) if none.
#'    - shared_operon: binary notation 0 for no overlap, 1 for overlap
#'
#' @return A data frame same as the df input, except gene pairs where either
#'    gene or both genes are not annotated to an operon are removed.
filter_missing_ko <- function(df) {
  keep <- lengths(df$ko_gene1) > 0 & lengths(df$ko_gene2) > 0
  df_filtered <- df[keep, ]
  message("Rows with missing KOs removed: ", sum(!keep))
  df_filtered
}

#' @description
#' Adds cumulative shared operon percentage column.
#' 
#' @param df A data frame where gene pairs are annotated to operons. Gene pairs 
#'    with shared operons are denoted, alongside a list of their shared operons. 
#'    Columns are described below:
#'    - gene1: gene name or locus tag of gene 1 of gene pair
#'    - gene2: gene name or locus tag of gene 2 of gene pair
#'    - correlation: similarity measure of gene pair
#'    - ko_gene1: operon ID or list of operon IDs annotated to gene 1. 
#'                NULL if none.
#'    - ko_gene2: operon ID or list of operon IDs annotated to gene 2
#'                NULL if none.
#'    - shared ko: operon ID or list of operon IDs annotated to both genes.
#'                character(0) if none.
#'    - shared_operon: binary notation 0 for no overlap, 1 for overlap
#'
#' @return A data frame same as the df input, with an additional column named
#'    cumulative_shared_pct. This column contains the cumulative percentage of
#'    shared operons up to that row. 
add_cumulative_pct <- function(df) {
  df$cumulative_shared_pct <- cumsum(df$shared_operon) / seq_along(df$shared_operon) * 100
  message("Cumulative percentage added")
  df
}