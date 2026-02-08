source("R/packages.R")
source("R/coexpression-functions.R")

# ================= LOAD DATA =================

genepairs <- readRDS("results/k12_genepairs.rds")
go_reference <- read_tsv("reference/derived/go_k12.tsv")

# ================= GO-ID MATCHES =================
filtered_pairs <- list(
  Original = filter_process_go(genepairs$Original,go_reference),
  MeanResiduals = filter_process_go(genepairs$MeanResiduals,go_reference),
  PC1Residuals = filter_process_go(genepairs$PC1Residuals,go_reference),
  CLR = filter_process_go(genepairs$CLR,go_reference),
  propr = filter_process_go(genepairs$propr,go_reference), 
  propr_z = filter_process_go(genepairs$propr_z,go_reference)
)

GO_mapped_pairs <- list(
  Original = get_process_go_ids(filtered_pairs$Original,go_reference,
                                cor_values = filtered_pairs$Original$correlation),
  MeanResiduals = get_process_go_ids(filtered_pairs$MeanResiduals,go_reference,
                                     cor_values = filtered_pairs$MeanResiduals$correlation),
  PC1Residuals = get_process_go_ids(filtered_pairs$PC1Residuals,go_reference,
                                    cor_values = filtered_pairs$PC1Residuals$correlation),
  CLR = get_process_go_ids(filtered_pairs$CLR,go_reference,
                           cor_values = filtered_pairs$CLR$correlation),
  propr = get_process_go_ids(filtered_pairs$propr,go_reference,
                             cor_values = filtered_pairs$propr$correlation), 
  propr_z = get_process_go_ids(filtered_pairs$propr_z,go_reference,
                               cor_values = filtered_pairs$propr_z$correlation)
)

genepairs_matches <- list(
  Original = compute_go_overlap(GO_mapped_pairs$Original),
  MeanResiduals = compute_go_overlap(GO_mapped_pairs$MeanResiduals),
  PC1Residuals = compute_go_overlap(GO_mapped_pairs$PC1Residuals),
  CLR = compute_go_overlap(GO_mapped_pairs$CLR),
  propr = compute_go_overlap(GO_mapped_pairs$propr), 
  propr_z = compute_go_overlap(GO_mapped_pairs$propr_z)
)

saveRDS(genepairs_matches,file="results/bacillus_GOmatches.rds")
