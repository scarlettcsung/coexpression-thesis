library(ggplot2)
library(ggridges)
library(dplyr)
library(tidyr)

reorder_matrix <- function(mat, sorted_names) {
  # 1. Identify genes that exist in both the sorted list and the matrix
  # This prevents errors if your sorted list has genes missing from the matrix
  common_genes <- sorted_names[sorted_names %in% rownames(mat)]
  
  if (length(common_genes) == 0) {
    stop("No matching gene names found between the sorted list and the matrix.")
  }
  
  # 2. Reorder rows and columns simultaneously
  return(mat[common_genes, common_genes])
}

# Usage:
# my_ordered_matrices <- reorder_matrices(my_similarity_matrices, my_sorted_list)

# ---------------------------------------------------------
# 1. Extraction Function
# ---------------------------------------------------------
extract_bin_correlations <- function(ordered_matrix_list) {
  
  results <- lapply(names(ordered_matrix_list), function(name) {
    mat <- ordered_matrix_list[[name]]
    n <- nrow(mat)
    dec <- floor(n / 10)
    
    # Indices
    low_idx <- 1:dec
    mid_idx <- floor((n - dec)/2 + 1):(floor((n - dec)/2 + 1) + dec - 1)
    top_idx <- (n - dec + 1):n
    
    # Create the list of vectors
    tmp_list <- list(
      "Bottom 10% to self"    = mat[low_idx, low_idx][upper.tri(mat[low_idx, low_idx])],
      "Middle 10% to self"    = mat[mid_idx, mid_idx][upper.tri(mat[mid_idx, mid_idx])],
      "Top 10% to self"       = mat[top_idx, top_idx][upper.tri(mat[top_idx, top_idx])],
      "Bottom 10% to Top 10%" = as.vector(mat[low_idx, top_idx])
    )
    
    # Use base R to stack and then dplyr to tidy up
    # This avoids the 'ind' not found issue
    df <- as.data.frame(utils::stack(tmp_list))
    
    # Manually set names to be safe
    colnames(df) <- c("correlation", "bin")
    df$method <- name
    return(df)
  })
  
  dplyr::bind_rows(results)
}

# ---------------------------------------------------------
# 2. Plotting Function
# ---------------------------------------------------------
plot_multidataset_collapsed <- function(combined_df, 
                                        method_order = NULL, 
                                        dataset_order = NULL) {
  
  # 1. Ordering
  if (!is.null(method_order)) {
    combined_df$method <- factor(combined_df$method, levels = method_order)
  }
  if (!is.null(dataset_order)) {
    combined_df$dataset <- factor(combined_df$dataset, levels = dataset_order)
  }
  
  # 2. Color Logic
  cor_colors <- c("#67001F","#B2182B","#D6604D","#F4A582","#FDDBC7","#FFFFFF",
                  "#D1E5F0","#92C5DE","#4393C3","#2166AC","#053061")
  color_fun <- grDevices::colorRamp(cor_colors)
  
  df_with_colors <- combined_df %>%
    group_by(dataset, method, bin) %>%
    mutate(grp_median = median(correlation, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
      fill_color = rgb(color_fun((grp_median + 1) / 2), maxColorValue = 255)
    )
  
  # 3. Plotting
  ggplot(df_with_colors, aes(x = correlation, y = method, fill = fill_color)) +
    geom_vline(xintercept = 0, linetype="dashed", color = "grey70") +
    
    geom_density_ridges(
      scale = 0.95, 
      quantile_lines = TRUE, 
      quantiles = 2, 
      color = "grey20",
      size = 0.2
    ) +
    
    scale_fill_identity() +
    
    # Grid: Datasets as Rows, Bins as Columns
    facet_grid(dataset ~ bin, switch = "y") + 
    
    coord_flip() +
    
    theme_ridges(center_axis_labels = TRUE) +
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
    labs(x = "Correlation", y = NULL) +
    
    theme(
      strip.background = element_rect(fill = "grey95"),
      strip.text = element_text(face = "bold", size = 10),
      panel.spacing = unit(0.2, "lines"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10)
    )
}

generate_median_table <- function(combined_df) {
  combined_df %>%
    group_by(dataset, method, bin) %>%
    summarise(median_cor = median(correlation, na.rm = TRUE), .groups = "drop") %>%
    # Pivot the bins to columns to make it easier to compare methods/datasets
    pivot_wider(names_from = bin, values_from = median_cor) %>%
    # Optional: Arrange for better readability
    arrange(dataset, method)
}
