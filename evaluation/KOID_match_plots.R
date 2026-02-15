source("R/packages.R")

# ================= DATA SETUP =================
ko_df <- readRDS("results/k12_komatches.rds")
n_genepairs <- floor(nrow(ko_df$Original) * 0.005)
total_mean <- ko_df$Original$cumulative_shared_pct[nrow(ko_df$Original)]

# BAR CHART
mean_at_subset <- sapply(names(ko_df),function(method) {
    df_method <- ko_df[[method]]
    return(df_method$cumulative_shared_pct[n_genepairs])
  })

mean_at_subset <- setNames(mean_at_subset, names(ko_df))

# Convert the named vector into a clean tibble
df_plot <- tibble(
  method = names(mean_at_subset),
  percent = mean_at_subset
) %>%
  mutate(
    diff_from_mean = percent - total_mean
  )

method_order <- c("Original","MeanResiduals","PC1Residuals",
                  "CLR","propr","propr_z")

# Set the factor levels (e.g., sorted by the percentage value)
df_plot <- df_plot %>%
  mutate(method = factor(method, levels = method_order))

# BINNED PLOT

binned_dataset_ko <- function(df,num_bins,dataset_mean) {
  comb_overlap <- do.call(cbind, lapply(df, function(x) x$shared_operon))
  colnames(comb_overlap) <- names(df)
  comb_overlap <- as.data.frame(comb_overlap)
  
  total_rows <- nrow(comb_overlap)
  N <- ceiling(total_rows / num_bins)
  comb_overlap$block <- ceiling(seq_len(total_rows) / N)
  comb_overlap$top_row <- pmin(comb_overlap$block * N, total_rows)
  
  summary_df <- comb_overlap %>%
    group_by(top_row) %>%
    summarize(
      Original = mean(Original, na.rm = TRUE)*100,
      MeanResiduals = mean(MeanResiduals, na.rm = TRUE)*100,
      PC1 = mean(PC1Residuals, na.rm = TRUE)*100,
      CLR = mean(CLR, na.rm = TRUE)*100,
      propr = mean(propr, na.rm = TRUE)*100,
      propr_z = mean(propr_z, na.rm = TRUE)*100
    )
  
  summary_df_long <- summary_df %>%
    tidyr::pivot_longer(cols = -top_row, names_to = "method", values_to = "percent")
  summary_df_long <- summary_df_long %>%
    mutate(
      diff_from_mean = (percent - dataset_mean)
    )
  return(summary_df_long)
}

num_bins = 1000
k12_summary <- binned_dataset_ko(ko_df,num_bins,total_mean)

# ================= PLOTTING =================

library(scales)
library(ggmagnify)

y_lab <- expression(Delta * " % GO Overlap")

# Bar Chart

barchart <- ggplot(df_plot, 
                   aes(x = method, y = diff_from_mean)) +
  geom_col(position = position_dodge(width = 0.9)) +
  theme_minimal() +
  labs(
    x = "Method",
    y = y_lab
  ) +
  scale_fill_brewer(palette = "Set1") +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

barchart

# Binned GO-ID % Matches

binned_k12 <- ggplot(k12_summary,
                     aes(x = top_row, y = diff_from_mean, color = method)) +
  geom_smooth(method = "loess", span = 0.01, se = FALSE, linewidth = .7) +
  labs(# title = "B. pseudomallei",
    x = "n gene pairs", 
    y = y_lab) +
  scale_x_continuous(labels = scales::label_comma(), # Added commas for readability
                     breaks = seq(0, max(k12_summary$top_row), by = 1000000)) +
  theme_minimal()+
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust = 1),
    # Move legend inside the plot
    legend.position = c(0.85, 0.65), 
    # Optional: Add a background/border to make it readable
    legend.background = element_rect(fill = "white", color = "lightgrey"))

# Combine plots
library(cowplot)

plot_grid(
  barchart, 
  binned_k12 + 
    geom_vline(xintercept = floor(nrow(ko_df$Original) * 0.005), 
               linetype = "dashed", color = "red", linewidth = 0.5) +
    geom_magnify(from = c(-0.1,500000,-0.1,5.8),
                 to = c(950000, 6500000, 1, 6)),
  labels = "AUTO",
  label_size = 14,
  label_fontface = "bold",
  nrow = 1,
  ncol = 2
)

# ================= PLOTTING (SUPPLEMENTARY) =================


ggplot(summary_df_long, aes(x = top_row, y = diff_from_mean)) +
  geom_line() +
  scale_x_continuous(
    breaks = seq(0, max_row, by = 1000000),
    labels = label_number(accuracy = 1) # show every 5th block
  ) +
  labs(x = "n gene pair (1000 bins) ",
       y = "Fraction of GO-matching pairs") +
  theme_minimal() +
  facet_wrap(~ method, nrow = 2)  # 2 rows, 3 plots per row

