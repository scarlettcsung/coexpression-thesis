source("R/packages.R")

dfs <- list(
  burk = readRDS("results/burkholderia_GOmatches.rds"),
  k12  = readRDS("results/k12_GOmatches.rds"),
  bac  = readRDS("results/bacillus_GOmatches.rds")
)

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

df_plot <- map2_dfr(mean_at_subset, names(mean_at_subset), ~{
  tibble(
    method = names(.x),
    percent = .x ,    # convert from fraction to percent
    dataset = .y
  )
})

df_plot <- df_plot %>%
  mutate(method = factor(method, levels = method_order))

ggplot(df_plot, aes(x = method, y = percent, fill = dataset)) +
  geom_col(position = position_dodge(width = 0.8)) +
  theme_minimal() +
  labs(
    title = "GO Overlap of Top 0.05% Correlated Gene Pairs",
    x = "Method",
    y = "Percent GO overlap",
    fill = "Dataset species"
  ) +
  scale_fill_brewer(palette = "Set2") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
