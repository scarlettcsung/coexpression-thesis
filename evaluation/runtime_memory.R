source("R/packages.R")
source("R/coexpression-functions.R")

library(bench)

path <- "data/k12_counts.csv" # File Path
raw_counts <- read.csv(path, header=TRUE, row.names = 1) # Load .csv
raw_counts <- raw_counts[order(colnames(raw_counts))]

# ================= CO-EXPRESSION ANALYSIS =================

# Filtered Counts
counts <- filterLowExpression(raw_counts) #Rows: Gene name, Columns: Conditions
counts_filtered <- counts[rowSums(counts) > 0, colSums(counts) > 0]
replicates = 2

bm1 <- bench::mark(
  Original       = corrs_Original(counts, replicates),
  MeanResiduals  = corrs_MeanResiduals(counts, replicates),
  PC1Residuals   = corrs_PC1Residuals(counts, replicates),
  iterations = 10,       # repeat each function 10 times
  check = FALSE,          # skip automatic equality checks
  gc = TRUE
)

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

bm <- bm_gc %>%
  as_tibble() %>%  # This is the magic line
  filter(!grepl("^gc$", as.character(expression)))

bm

bm %>%
  dplyr::select(expression, time) %>%
  group_by(expression) %>%
  summarise(
    mean_time = mean(as.numeric(as_bench_time(unlist(time)), units = "ms")),
    sd_time   = sd(as.numeric(as_bench_time(unlist(time)), units = "ms"))
  )


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

library(patchwork)
combined_plot <- p1 + p2 + 
  plot_layout(widths = c(1.2, 1))

# Display
combined_plot
