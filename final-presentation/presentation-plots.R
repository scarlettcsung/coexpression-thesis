source("R/packages.R")
source("R/plotting-functions.R")

# This code is a mess but it's just plots for the presentation

# File path
path <- "data/BURK_counts.csv"

# Load .csv
raw_counts <- read.csv(path, header=TRUE, row.names = 1) 
# Orders samples by alphabetical order
raw_counts <- raw_counts[order(colnames(raw_counts))] 

burk_sorted_names <- readRDS("reference/derived/burkholderia_genenames_sorted.rds")
# ================= CO-EXPRESSION ANALYSIS =================

# Filtering Counts -------------------------------------------------------------
counts <- filterLowExpression(raw_counts) 

tmm <- logTMM(counts,replicates = 3)
sorted <- tmm[order(rowMeans(tmm)), ]

# Just PCC ---------------------------------------------------------------------
pcc_og <- WGCNA::cor(t(sorted))
pcc_df <- data.frame(x=pcc_og[upper.tri(pcc_og)])

# Compute stats separately
median_x <- median(pcc_df$x, na.rm = TRUE)
q1 <- quantile(pcc_df$x, 0.25, na.rm = TRUE)
q3 <- quantile(pcc_df$x, 0.75, na.rm = TRUE)

# Plot Full Plot ---------------------------------------------------------------
ggplot(pcc_df, aes(x = x)) + 
  geom_density(fill = "white", color = "black", alpha = 0.7) +
  geom_vline(xintercept = median_x, color = "red", linetype = "dashed", size = 1) +
  geom_vline(xintercept = q1, color = "blue", linetype = "dotted", size = 1) +
  geom_vline(xintercept = q3, color = "blue", linetype = "dotted", size = 1) + 
  labs(title = "Distribution of Correlation Coefficients\n",
         x = "Correlation", y = "Density") +
  theme_minimal()

reordered <- reorder_matrix(pcc_og, burk_sorted_names)

extract_single_bin_correlations <- function(mat) {
  
  n <- nrow(mat)
  dec <- floor(n / 10)
  
  # Indices
  low_idx <- 1:dec
  mid_idx <- floor((n - dec)/2 + 1):(floor((n - dec)/2 + 1) + dec - 1)
  top_idx <- (n - dec + 1):n
  
  # Extract correlations for each bin
  bin_list <- list(
    "Bottom 10% to self"    = mat[low_idx, low_idx][upper.tri(mat[low_idx, low_idx])],
    "Middle 10% to self"    = mat[mid_idx, mid_idx][upper.tri(mat[mid_idx, mid_idx])],
    "Top 10% to self"       = mat[top_idx, top_idx][upper.tri(mat[top_idx, top_idx])],
    "Bottom 10% to Top 10%" = as.vector(mat[low_idx, top_idx])
  )
  
  # Stack into a tidy data frame
  df <- as.data.frame(stack(bin_list))
  colnames(df) <- c("correlation", "bin")
  
  return(df)
}


binned <- extract_single_bin_correlations(reordered)

plot_single_dataset <- function(df) {
  
  # Compute median per bin (optional, for coloring)
  cor_colors <- c("#67001F","#B2182B","#D6604D","#F4A582","#FDDBC7","#FFFFFF",
                  "#D1E5F0","#92C5DE","#4393C3","#2166AC","#053061")
  color_fun <- grDevices::colorRamp(cor_colors)
  
  df <- df %>%
    group_by(bin) %>%
    mutate(grp_median = median(correlation, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(
      fill_color = rgb(color_fun((grp_median + 1) / 2), maxColorValue = 255)
    )
  
  # Plot
  ggplot(df, aes(x = correlation, y = bin, fill = fill_color)) +
    geom_vline(xintercept = 0, linetype="dashed", color = "grey70") +
    geom_density_ridges(
      scale = 0.95, 
      quantile_lines = TRUE, 
      quantiles = 2, 
      color = "grey20",
      size = 0.2
    ) +
    scale_fill_identity() +
    coord_flip() +
    scale_y_discrete(expand = c(0, 0)) +   # <- remove the left/right gap
    theme_ridges(center_axis_labels = TRUE) +
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
    labs(x = "Correlation", y = "Density") +
    theme(
      axis.text.x = element_text(hjust = -.2, size = 10),
      axis.text.y = element_text(size = 10)
    )
}

plot_single_dataset(binned)

# Re-order matrices 
matrices_burk <- readRDS("results_ignore/burkholderia_matrices.rds")
burk_sorted_names <- readRDS("reference/derived/burkholderia_genenames_sorted.rds")
burk_ordered_matrices <- lapply(matrices_burk, function(m) {
  reorder_matrix(m, burk_sorted_names)
})

spqn_reordered <- reorder_matrix(matrices_burk$Original, burk_sorted_names)
spqn_binned <- extract_single_bin_correlations(spqn_reordered)
plot_single_dataset(spqn_binned)

mr_reordered <- reorder_matrix(matrices_burk$MeanResiduals, burk_sorted_names)
mr_binned <- extract_single_bin_correlations(mr_reordered)
plot_single_dataset(mr_binned)

pc1_reordered <- reorder_matrix(matrices_burk$PC1Residuals, burk_sorted_names)
pc1_binned <- extract_single_bin_correlations(pc1_reordered)
plot_single_dataset(pc1_binned)
