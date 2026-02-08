source("R/packages.R")
source("R/coexpression-functions.R")
# ================= DATA SETUP =================
matrices <- readRDS("results/k12_matrices.rds")

# ================= HIERARCHAL CLUSTERING OVERVIEW =================

distances <- lapply(matrices, function(mat) 1 - mat)

library(dynamicTreeCut)
geneTrees <- lapply(distances, function(diss) {
  hclust(as.dist(diss), method = "average")
})

library(dendextend)
library(colorspace)

cluster_results <- lapply(names(geneTrees), function(nm) {
  
  tree <- geneTrees[[nm]]
  
  # Run dynamic tree cut
  clusters <- cutreeDynamic(
    dendro = tree,
    minClusterSize = 30,   # minimum genes per cluster
    deepSplit = 1,
    method = "tree"
  )
  
  names(clusters) <- tree$labels
  
  # Optional: drop unassigned clusters (0) if needed
  clusters_filtered <- clusters
  clusters_filtered[clusters_filtered == 0] <- NA
  
  # Count clusters after filtering
  n_clusters <- length(unique(na.omit(clusters_filtered)))
  
  
  # Return results
  list(
    clusters = clusters_filtered,
    n_clusters = n_clusters,
    cluster_sizes = table(clusters_filtered)
  )
})

par(mfrow = c(2, 3))  # adjust to number of methods

for (nm in names(geneTrees)) {
  
  tree <- geneTrees[[nm]]
  
  # Dynamic tree cut
  clusters <- cutreeDynamic(tree, 
                            minClusterSize = 30, 
                            deepSplit = 1, 
                            method = "tree")
  clusters[clusters == 0] <- NA
  
  # Make cluster numbers contiguous
  clusters_noNA <- clusters[!is.na(clusters)]
  clusters_contig <- as.numeric(factor(clusters_noNA))
  clusters_colored <- rep(NA, length(clusters))
  clusters_colored[!is.na(clusters)] <- clusters_contig
  
  # Convert to dendrogram
  dend <- as.dendrogram(tree)
  dend <- reorder(dend, wts = clusters_colored)
  dend <- color_branches(dend, clusters = clusters_colored)
  
  plot(dend, main = nm, ylab = "1 - cor", xlab = "", leaflab = "none")
}

names(cluster_results) <- names(geneTrees)
data.frame(
  n_clusters = sapply(cluster_results, `[[`, "n_clusters")
)
# ================= FILTER CLUSTERS =================

for (nm in names(cluster_results)) {
  cat("\nMethod:", nm)
  print(table(cluster_results[[nm]]$clusters))
}

# IGNORE FOR NOW
min_size <- 30
filtered_clusters <- lapply(cluster_results, function(res) {
  cl <- res$clusters
  tab <- table(cl)
  
  keep <- names(tab)[tab >= min_size]
  cl[cl %in% keep]
})

print(k)
sapply(filtered_clusters, function(cl) {
  length(unique(cl))/k
})
# STOP IGNORING

saveRDS(cluster_results,"evaluation/k12_clusters.rds")

# ================= GO ENRICHMENT =================

library(dplyr)
library(tidyr)
library(stringr)

# Make reference database to match locus tag to gene name for GO Enrichment
gnr <- read_tsv("reference/derived/k12_gene_name_ref.tsv")
gnr_long <- gnr %>%
  dplyr::select(`Gene Names (ordered locus)`, `Gene Names (primary)`) %>%
  tidyr::separate_rows(`Gene Names (ordered locus)`, sep = " ") %>%
  dplyr::rename(locus_tag = `Gene Names (ordered locus)` ,
                symbol = `Gene Names (primary)`) %>%
  dplyr::filter(locus_tag != "")   # drop empty strings

clusters_symbols <- lapply(names(cluster_results), function(nm) {
  
  cluster_vec <- cluster_results[[nm]]$clusters  # vector of cluster assignments
  names(cluster_vec) <- geneTrees[[nm]]$labels
  
  cluster_df <- data.frame(
    Gene = names(cluster_vec),
    Cluster = cluster_vec,
    Method = nm
  )
  
  # Join with your gene info
  cluster_df <- cluster_df %>%
    left_join(gnr_long, by = c("Gene" = "locus_tag"))
  
  cluster_df
})

names(clusters_symbols) <- names(cluster_results)

all_genes <- data.frame(
  Gene = geneTrees$Original$labels   
) %>%
  left_join(gnr_long, by = c("Gene" = "locus_tag")) %>%
  pull(symbol)

library(clusterProfiler)
library(org.EcK12.eg.db)
library(AnnotationDbi) 
library(GOSemSim)
library(enrichplot)

ego_all <- lapply(names(clusters_symbols), function(nm) {
  cl_df <- clusters_symbols[[nm]]
  cl_list <- split(cl_df$symbol, cl_df$Cluster)
  
  compareCluster(
    geneCluster = cl_list,
    universe = all_genes,
    fun = "enrichGO",
    OrgDb = org.EcK12.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    pAdjustMethod = "BH",
    qvalueCutoff = 0.01
  )
})

names(ego_all) <- names(clusters_symbols)

saveRDS(ego_all,"evaluation/k12_ego.rds")
ego_all <- readRDS("evaluation/k12_ego.rds")

ego_simpl <- lapply(names(ego_all),function(nm) {
  simplified <- simplify(ego_all[[nm]],
                         cutoff=0.07,
                         by='p.adjust',
                         measure='Wang')
})
names(ego_simpl) <- names(ego_all)
  
# ego_simpl <- pairwise_termsim(ego_simpl)
# emapplot(ego_simpl)

ego_df <- lapply(ego_simpl,function(x) as.data.frame(x))

saveRDS(ego_simpl,"evaluation/k12_ego_simpl.rds")
ego_df <- readRDS("evaluation/k12_ego_simpl.rds")
ego_df <- as.data.frame(ego_df)

view(ego_df$Original)
view(ego_df$MeanResiduals)
view(ego_df$PC1Residuals)
view(ego_df$CLR)
view(ego_df$propr)
view(ego_df$propr_z)

ego_df <- lapply(ego_df,function(x) as.data.frame(x))

ego_top <- lapply(ego_df, function(df) {
  df %>%
    group_by(Cluster) %>%
    filter(GeneRatio == max(GeneRatio, na.rm = TRUE)) %>%
    ungroup()
})

viewer <- ego_top$Original
view(viewer[,c("Cluster","GeneRatio")])
view(viewer[,c("Cluster","ID","Description","GeneRatio","geneID")])

go_lists <- lapply(ego_top, function(df) df$ID)
names(go_lists) <- names(ego_top)
# All unique GO terms
all_go <- unique(unlist(go_lists))

# Rows = GO terms, columns = methods, TRUE if present
go_matrix <- sapply(go_lists, function(go) all_go %in% go)
rownames(go_matrix) <- all_go
go_matrix <- as.data.frame(go_matrix)
go_matrix$method_count <- rowSums(go_matrix)
go_matrix$methods <- apply(go_matrix[, names(go_lists)], 1, function(x) {
  paste(names(x)[x], collapse = ", ")
})

# Create a mapping of ID → Description
# Take the first description from any method that contains it
go_desc <- lapply(ego_df, function(df) {
  df %>% select(ID, Description)
}) %>%
  bind_rows() %>%
  distinct(ID, .keep_all = TRUE)

# Join with go_matrix
go_matrix <- go_matrix %>%
  rownames_to_column(var = "ID") %>%
  left_join(go_desc, by = "ID") %>%
  column_to_rownames(var = "ID")

view(go_matrix)

method <- "propr_z"
cluster_num <- 1
clusters <- cluster_results[[method]]$clusters

genes_to_plot <-names(clusters)[clusters == cluster_num]
genes_to_plot <- genes_to_plot[!is.na(genes_to_plot)]
mat_subset <- matrices[[method]][genes_to_plot,genes_to_plot]

library(pheatmap)

par(mfrow = c(1, 3))  # adjust to number of methods

pheatmap(mat_subset,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = paste("Cluster", cluster_num, "Heatmap -", method)
)

long_df <- go_matrix %>%
  select(Description, methods) %>%
  separate_rows(methods, sep = ",\\s*")

wide_df <- long_df %>%
  mutate(value = 1L) %>%
  pivot_wider(
    names_from = methods,
    values_from = value,
    values_fill = 0
  ) %>%
  rowwise() %>%
  mutate(method_sum = sum(c_across(where(is.numeric)))) %>%
  ungroup()