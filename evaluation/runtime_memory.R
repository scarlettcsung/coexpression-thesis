source("R/packages.R")
source("R/coexpression-functions.R")

librarian::shelf(
  bench,patchwork
)

path <- "data/k12_counts.csv" # File Path
raw_counts <- read.csv(path, header=TRUE, row.names = 1) # Load .csv
raw_counts <- raw_counts[order(colnames(raw_counts))]

# ================= CO-EXPRESSION ANALYSIS =================

# Filtered Counts
counts <- filterLowExpression(raw_counts) #Rows: Gene name, Columns: Conditions
replicates = 2

# gc to make results more consistent
# bms are broken in 2 as it can be a bit heavy for memory

gc()

bm1 <- bench::mark(
  Original       = corrs_Original(counts, replicates),
  MeanResiduals  = corrs_MeanResiduals(counts, replicates),
  PC1Residuals   = corrs_PC1Residuals(counts, replicates),
  iterations = 10,       # repeat each function 10 times
  check = FALSE,          # skip automatic equality checks
  gc = TRUE
)

gc()

bm2 <- bench::mark(
  CLR            = corrs_CLR(counts),
  propr          = get_propr(counts),
  propr_z        = get_zpropr(counts),
  iterations = 10,       # repeat each function 10 times
  check = FALSE,          # skip automatic equality checks
  gc = TRUE
)

bm_gc <- bind_rows(bm1,bm2)

bm_gc

# combine
bm <- bm_gc %>%
  as_tibble() %>%  # This is the magic line
  filter(!grepl("^gc$", as.character(expression)))

bm

# view mean and sds
bm %>%
  dplyr::select(expression, time) %>%
  group_by(expression) %>%
  summarise(
    mean_time = mean(as.numeric(as_bench_time(unlist(time)), units = "ms")),
    sd_time   = sd(as.numeric(as_bench_time(unlist(time)), units = "ms"))
  )

saveRDS(bm,"evaluation/eval_results/comp_benchmark_results.rds")

# Plotting
plot_data <- bm %>%
  dplyr::select(expression,time) %>%
  unnest(time) %>%
  mutate(
    expression = as.character(expression),
    time_sec = as.numeric(time)
  )

method_order <- c("Original", "MeanResiduals", "PC1Residuals", 
                  "CLR", "propr", "propr_z")

plot_data <- plot_data %>%
  mutate(expression = factor(expression, levels = my_order))

p1 <- ggplot(plot_data, aes(x = fct_rev(expression), y = time_sec)) +
  geom_boxplot() +
  coord_flip()+
  labs(x = "Method", y = "Runtime (seconds)")

bm_for_plot <- bm %>%
  as_tibble() %>%
  mutate(
    median_num = as.numeric(median),
    mem_gb = as.numeric(mem_alloc) / 1024^3,
    # CRITICAL: Don't use as.character(), use factor() with your levels
    expression = factor(expression, levels = method_order)
  )

p2 <- ggplot(bm_for_plot, aes(x = fct_rev(expression), y = mem_gb)) +
  geom_col(fill = "slategray3", width = 0.7) +
  coord_flip() +
  theme_minimal() +
  labs(x = NULL, y = "Memory (GB)") +
  theme(axis.text.y = element_blank())


combined_plot <- p1 + p2 + 
  plot_layout(widths = c(1.2, 1))

# Display
combined_plot

bm <- readRDS("evaluation/eval_results/comp_benchmark_results.rds")

# Friedman test

bm_friedman <- bm %>%
  dplyr::select(expression, time) %>%
  tidyr::unnest(time) %>%                  # expand the list-column of times
  mutate(
    time_s = as.numeric(time)
    ) %>%
  dplyr::select(expression, time_s)

friedman_mat_df <- bm_friedman %>%
  group_by(expression) %>%
  mutate(iter = row_number()) %>%
  tidyr::pivot_wider(names_from = expression, values_from = time_s)

friedman_mat <- as.matrix(friedman_mat_df[,-1])
friedman <- friedman.test(friedman_mat)  

n_iter <- nrow(friedman_mat)
n_method <- ncol(friedman_mat)
x <- as.vector(friedman_mat)
g <- rep(colnames(friedman_mat), each = n_iter) 
b <- rep(1:n_iter, times = n_method)  
xDF <- data.frame(x = x, g = factor(g), b = b)

frdAllPairsNemenyiTest(xDF$x[o], groups = xDF$g[o], blocks = xDF$b[o])
