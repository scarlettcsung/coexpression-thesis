source("R/packages.R")

# ================= DATA SETUP =================
dfs <- list(
  B_pseudomallei = readRDS("results/burkholderia_GOmatches.rds"),
  E_coli  = readRDS("results/k12_GOmatches.rds"),
  B_subtilis  = readRDS("results/bacillus_GOmatches.rds")
)

# BAR CHART
mean_at_subset <- lapply(names(dfs),function(nm) {
  df <- dfs[[nm]]
  n_genepairs <- floor(nrow(df$Original) * 0.005)
  
  sapply(names(df),function(method) {
    df_method <- df[[method]]
    df_method$GO_match_cumulative[n_genepairs]
  })
})
mean_at_subset <- setNames(mean_at_subset, names(dfs))

method_order <- names(mean_at_subset[[1]])
dataset_order <- c("B_pseudomallei","E_coli","B_subtilis")

df_plot <- map2_dfr(mean_at_subset, names(mean_at_subset), ~{
  tibble(
    method = names(.x),
    percent = .x ,    # convert from fraction to percent
    dataset = .y
  )
})

full_means <- c("B_pseudomallei" = 2.87786, 
                "E_coli" = 4.87045, 
                "B_subtilis" = 5.65051)

df_plot <- df_plot %>%
  mutate(
    dataset_mean = full_means[dataset],
    diff_from_mean = percent - dataset_mean
  )

df_plot <- df_plot %>%
  mutate(method = factor(method, levels = method_order),
         dataset = factor(dataset, levels = dataset_order))

# BINNED PLOT

binned_dataset <- function(df,num_bins,dataset_mean) {
  comb_overlap <- do.call(cbind, lapply(df, function(x) x$GO_overlap))
  colnames(comb_overlap) <- names(df)
  comb_overlap <- as.data.frame(comb_overlap)
  
  total_rows <- nrow(comb_overlap)
  N <- ceiling(total_rows / num_bins)
  comb_overlap$block <- ceiling(seq_len(total_rows) / N)
  comb_overlap$top_row <- pmin(comb_overlap$block * N, total_rows)
  
  summary_df <- comb_overlap %>%
    group_by(top_row) %>%
    summarize(
      Original = mean(Original, na.rm = TRUE),
      MeanResiduals = mean(MeanResiduals, na.rm = TRUE),
      PC1 = mean(PC1Residuals, na.rm = TRUE),
      CLR = mean(CLR, na.rm = TRUE),
      propr = mean(propr, na.rm = TRUE),
      propr_z = mean(propr_z, na.rm = TRUE)
    )

  summary_df_long <- summary_df %>%
    tidyr::pivot_longer(cols = -top_row, names_to = "method", values_to = "percent")
  summary_df_long <- summary_df_long %>%
    mutate(
      diff_from_mean = (percent - dataset_mean)*100
  )
  return(summary_df_long)
}

num_bins = 1000

burk_summary <- binned_dataset(dfs$B_pseudomallei,num_bins,0.0287786)
k12_summary <- binned_dataset(dfs$E_coli,num_bins,0.0487045)
bac_summary <- binned_dataset(dfs$B_subtilis,num_bins,0.0565051)

# ================= PLOTTING =================

library(scales)
library(ggmagnify)

y_lab <- expression(Delta * " % GO Overlap")

# Bar Chart

barchart <- ggplot(df_plot, 
                   aes(x = method, y = diff_from_mean, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.9)) +
  theme_minimal() +
  labs(
    x = "Method",
    y = y_lab,
    fill = "Dataset species"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

# Binned GO-ID % Matches

binned_burk <- ggplot(burk_summary, 
                      aes(x = top_row, y = diff_from_mean, color = method)) +
  geom_smooth(method = "loess", span = 0.01, se = FALSE, linewidth = .7) +
  labs(# title = "B. pseudomallei",
       x = "n gene pairs", 
       y = y_lab) +
  scale_x_continuous(labels = scales::label_comma(), # Added commas for readability
                     breaks = seq(0, max(burk_summary$top_row), by = 500000)) +
  scale_y_continuous(labels = scales::label_comma(), # Added commas for readability
                     breaks = seq(0, 22, by = 5)) +
  theme_minimal()+
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust = 1),
    # Move legend inside the plot
    legend.position = c(0.85, 0.65), 
    # Optional: Add a background/border to make it readable
    legend.background = element_rect(fill = "white", color = "lightgrey"))


binned_k12 <- ggplot(k12_summary, 
                      aes(x = top_row, y = diff_from_mean, color = method)) +
  geom_smooth(method = "loess", span = 0.01, se = FALSE, linewidth = .7) +
  labs(# title = "E. coli",
       x = "n gene pairs", 
       y = y_lab) +
  scale_x_continuous(labels = scales::label_comma(), # Added commas for readability
                     breaks = seq(0, max(k12_summary$top_row), by = 1000000))+
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none")

binned_bac <- ggplot(bac_summary, 
                      aes(x = top_row, y = diff_from_mean, color = method)) +
  geom_smooth(method = "loess", span = 0.01, se = FALSE, linewidth = .7) +
  labs(# title = "B. subtilis",
       x = "n gene pairs", 
       y = y_lab) +
  scale_x_continuous(labels = scales::label_comma(), # Added commas for readability
                     breaks = seq(0, max(bac_summary$top_row), by = 500000)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust = 1), 
    legend.position = "none")

# Combine plots
library(cowplot)

plot_grid(
  barchart, 
  binned_burk + 
    geom_vline(xintercept = floor(nrow(dfs$B_pseudomallei$Original) * 0.005), 
               linetype = "dashed", color = "red", linewidth = 0.5) +
    geom_magnify(from = c(-0.2,300000,-1,22),
                 to = c(420000, 2000000, 4, 24)),
  binned_k12 + 
    geom_vline(xintercept = floor(nrow(dfs$E_coli$Original) * 0.005), 
               linetype = "dashed", color = "red", linewidth = 0.5) +
    geom_magnify(from = c(-0.1,300000,-0.5,17),
                 to = c(750000, 4200000, 4, 18)),
  binned_bac + 
    geom_vline(xintercept = floor(nrow(dfs$B_subtilis$Original) * 0.005), 
               linetype = "dashed", color = "red", linewidth = 0.5) +
    geom_magnify(from = c(-0.1,200000,-1,36),
                 to = c(450000, 2250000, 8, 40)),
  labels = "AUTO",
  label_size = 14,
  label_fontface = "bold",
  nrow = 2,
  ncol = 2
)

# ================= PLOTTING (ADDITIONAL) =================

combined_burk <- do.call(cbind, lapply(dfs$B_pseudomallei, function(x) x$GO_match_cumulative))
combined_k12 <- do.call(cbind, lapply(df_list, function(x) x$GO_match_cumulative))
combined_bac <- do.call(cbind, lapply(df_list, function(x) x$GO_match_cumulative))

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

