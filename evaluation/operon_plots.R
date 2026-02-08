source("R/packages.R")
source("R/coexpression-functions.R")

ko_matches <- saveRDS(ko_matches,file="results/k12_komatches.rds")

plot_df <- do.call(
  rbind,
  lapply(names(ko_matches), function(m) {
    df <- ko_matches[[m]]
    data.frame(
      method = m,
      pair_index = seq_len(nrow(df)),
      cumulative_shared_pct = df$cumulative_shared_pct
    )
  })
)

library(ggplot2)

ggplot(plot_df,
       aes(x = pair_index,
           y = cumulative_shared_pct,
           color = method)) +
  geom_line(linewidth = 1) +
  labs(
    x = "Gene pair rank",
    y = "Cumulative % of shared operons",
    color = "Method"
  ) +
  theme_minimal()

max_n <- nrow(ko_matches$Original)
max_n <- 50000

plot(
  x = c(1, max_n),
  y = c(0, 100),
  type = "n",
  xlab = "Gene pair rank",
  ylab = "Cumulative % of shared operons"
)

i <- 1

for (m in names(ko_matches)) {
  df <- ko_matches[[m]]
  
  lines(
    x = seq_len(nrow(df)),
    y = df$cumulative_shared_pct,
    col = cols[i],
    lwd = 2
  )
  
  i <- i + 1
}

legend(
  "topright",
  legend = names(ko_matches),
  col = cols,
  lwd = 2,
  cex = 0.8
)
