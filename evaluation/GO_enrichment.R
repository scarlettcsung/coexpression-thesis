source("R/packages.R")

# ================= LOAD DEPENDENCIES =================
# Install and load other required packages
librarian::shelf(
  tidyr, stringr,
  dynamicTreeCut,
  dendextend, colorspace
)

# Install if needed
BiocManager::install(c("clusterProfiler","AnnotationDbi",
                       "org.EcK12.eg.db","enrichplot",
                       "GO.db"))

# ================= DATA SETUP =================
# Load matrices
matrices <- readRDS("results_ignore/k12_matrices.rds")

# ================= HIERARCHAL CLUSTERING OVERVIEW =================
# Distance measure for histogram
distances <- lapply(matrices, function(mat) 1 - mat)

# Make histogram
geneTrees <- lapply(distances, function(diss) {
  hclust(as.dist(diss), method = "average")
})

# Clustering with cutreeDynamic
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

# Clustering with cutreeDynamic and plotting
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

# Returning method names
names(cluster_results) <- names(geneTrees)
data.frame(
  n_clusters = sapply(cluster_results, `[[`, "n_clusters")
)

for (nm in names(cluster_results)) {
  cat("\nMethod:", nm)
  print(table(cluster_results[[nm]]$clusters))
}
saveRDS(cluster_results,"evaluation/k12_clusters.rds")

# ================= GO ENRICHMENT =================

# Make reference database to match locus tag to gene name for GO Enrichment
gnr <- read_tsv("reference/derived/k12_gene_name_ref.tsv")
gnr_long <- gnr %>%
  dplyr::select(`Gene Names (ordered locus)`, `Gene Names (primary)`) %>%
  tidyr::separate_rows(`Gene Names (ordered locus)`, sep = " ") %>%
  dplyr::rename(locus_tag = `Gene Names (ordered locus)` ,
                symbol = `Gene Names (primary)`) %>%
  dplyr::filter(locus_tag != "")   # drop empty strings

# Changing to names that are consistent with GO database
# nm is method names!
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

# Making "Universe" data frame
all_genes <- data.frame(
  Gene = geneTrees$Original$labels   
) %>%
  left_join(gnr_long, by = c("Gene" = "locus_tag")) %>%
  pull(symbol)

# Run GO-enrichment for each method
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

# Simplify terms
ego_simpl <- lapply(names(ego_all),function(nm) {
  simplified <- simplify(ego_all[[nm]],
                         cutoff=0.07,
                         by='p.adjust',
                         measure='Wang')
})
names(ego_simpl) <- names(ego_all)

# Make into dfs
ego_df <- lapply(ego_simpl,function(x) as.data.frame(x))

saveRDS(ego_simpl,"evaluation/k12_ego_simpl.rds")
ego_df <- readRDS("evaluation/k12_ego_simpl.rds")
ego_df <- as.data.frame(ego_df)

ego_df <- lapply(ego_df,function(x) as.data.frame(x))

# Top terms as representative terms
ego_top <- lapply(ego_df, function(df) {
  df %>%
    group_by(Cluster) %>%
    filter(GeneRatio == max(GeneRatio, na.rm = TRUE)) %>%
    ungroup()
})

# To view clusters
viewer <- ego_top$Original
view(viewer[,c("Cluster","GeneRatio")])
view(viewer[,c("Cluster","ID","Description","GeneRatio","geneID")])

# Organizing as dataframes

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
