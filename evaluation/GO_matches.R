source("R/packages.R")
source("R/GO-overlap-functions.R")

# Usage: run in IDE, this script matches genes to GO-IDs (process) and 
#        calculates GO-ID overlap. 

# ================= LOAD DATA =================

# Gene pairs were not uploaded to GitHub repository due to file size. Generate
# this with 01-coexpression-analysis.R.
genepairs <- readRDS("results_ignore/k12_genepairs.rds")
go_reference <- read_tsv("reference/derived/go_k12.tsv")

# ================= GO-ID MATCHES =================
# Filtering gene pairs of all methods
# Both genes of gene pairs must match to at least 1 GO-ID
filtered_pairs <- list(
  Original = filter_process_go(genepairs$Original,go_reference),
  MeanResiduals = filter_process_go(genepairs$MeanResiduals,go_reference),
  PC1Residuals = filter_process_go(genepairs$PC1Residuals,go_reference),
  CLR = filter_process_go(genepairs$CLR,go_reference),
  propr = filter_process_go(genepairs$propr,go_reference), 
  propr_z = filter_process_go(genepairs$propr_z,go_reference)
)

# Maps GO-IDs to genes
GO_mapped_pairs <- list(
  Original = get_process_go_ids(
    filtered_pairs$Original,go_reference,
    cor_values = filtered_pairs$Original$correlation),
  MeanResiduals = get_process_go_ids(
    filtered_pairs$MeanResiduals,go_reference,
    cor_values = filtered_pairs$MeanResiduals$correlation),
  PC1Residuals = get_process_go_ids(
    filtered_pairs$PC1Residuals,go_reference,
    cor_values = filtered_pairs$PC1Residuals$correlation),
  CLR = get_process_go_ids(
    filtered_pairs$CLR,go_reference,
    cor_values = filtered_pairs$CLR$correlation),
  propr = get_process_go_ids(
    filtered_pairs$propr,go_reference,
    cor_values = filtered_pairs$propr$correlation), 
  propr_z = get_process_go_ids(
    filtered_pairs$propr_z,go_reference,
    cor_values = filtered_pairs$propr_z$correlation)
)

# Computes whether gene pairs have GO-overlap, and the cumulative overlaps
genepairs_matches <- list(
  Original = compute_go_overlap(GO_mapped_pairs$Original),
  MeanResiduals = compute_go_overlap(GO_mapped_pairs$MeanResiduals),
  PC1Residuals = compute_go_overlap(GO_mapped_pairs$PC1Residuals),
  CLR = compute_go_overlap(GO_mapped_pairs$CLR),
  propr = compute_go_overlap(GO_mapped_pairs$propr), 
  propr_z = compute_go_overlap(GO_mapped_pairs$propr_z)
)

# Saves GO matches to .rds file for future use. These files are too big to be 
# included in GitHub repository but can be generated here.
saveRDS(genepairs_matches,file="results_ignore/k12_GOmatches.rds")
