source("R/packages.R")
source("R/coexpression-functions.R")

# ================= LOAD DATA =================

genepairs <- readRDS("results/k12_genepairs.rds")
ko_reference <- read_csv("reference/derived/k12_operon_reference.csv")

# ================= GO-ID MATCHES =================
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

filter_missing_ko <- function(df) {
  keep <- lengths(df$ko_gene1) > 0 & lengths(df$ko_gene2) > 0
  df_filtered <- df[keep, ]
  message("Rows with missing KOs removed: ", sum(!keep))
  df_filtered
}

add_cumulative_pct <- function(df) {
  df$cumulative_shared_pct <- cumsum(df$shared_operon) / seq_along(df$shared_operon) * 100
  message("Cumulative percentage added")
  df
}

shared_ko <- lapply(genepairs, compute_shared_ko,
                             gene_to_ko = ko_reference)

saveRDS(shared_ko,file="results/k12_shared.rds")

filtered_shared_ko <- lapply(shared_ko, filter_missing_ko)

ko_matches <- lapply(filtered_shared_ko, add_cumulative_pct)

saveRDS(ko_matches,file="results/k12_komatches.rds")
