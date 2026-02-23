source("R/packages.R")
librarian::shelf(RColorBrewer)

# ================= DATA SETUP =================
dfs <- list(
  B_pseudomallei = readRDS("results/burkholderia_GOmatches.rds"),
  E_coli  = readRDS("results/k12_GOmatches.rds"),
  B_subtilis  = readRDS("results/bacillus_GOmatches.rds")
)

# ===================== Setup =====================
step <- 10  # thinning step

# Matrices
comb_list <- list(
  B_pseudomallei = do.call(cbind, lapply(dfs$B_pseudomallei, function(x) x$GO_match_cumulative)),
  E_coli         = do.call(cbind, lapply(dfs$E_coli, function(x) x$GO_match_cumulative)),
  B_subtilis     = do.call(cbind, lapply(dfs$B_subtilis, function(x) x$GO_match_cumulative))
)

corr_list <- list(
  B_pseudomallei = do.call(cbind, lapply(dfs$B_pseudomallei, function(x) x$correlation)),
  E_coli         = do.call(cbind, lapply(dfs$E_coli, function(x) x$correlation)),
  B_subtilis     = do.call(cbind, lapply(dfs$B_subtilis, function(x) x$correlation))
)

# Colors and method names
cols <- brewer.pal(6,"Dark2") 
methods <- colnames(comb_list[[1]])

# Desired x-axis ranges per dataset
x_ranges <- list(
  B_pseudomallei = c(1, 100000),
  E_coli         = c(1, 100000),
  B_subtilis     = c(1, 100000)
)

# ===================== Subset & thin =====================
x_list <- list()
comb_sub_list <- list()
corr_sub_list <- list()

for (name in names(comb_list)) {
  mat <- comb_list[[name]]
  corr <- corr_list[[name]]
  
  # Original x indices (gene pair numbers)
  x_full <- 1:nrow(mat)
  
  # Thinning
  thin_idx <- seq(1, length(x_full), by = step)
  
  # Subset to desired x-range
  thin_idx <- thin_idx[x_full[thin_idx] <= x_ranges[[name]][2]]
  
  # Subset both x and matrices using same indices
  x_list[[name]] <- x_full[thin_idx]
  comb_sub_list[[name]] <- mat[thin_idx, , drop = FALSE]
  corr_sub_list[[name]] <- corr[thin_idx, , drop = FALSE]
}

# ===================== Plot =====================
par(mfrow = c(2,3), mar = c(4,4,2,1),bg="white")

# Top row: cumulative % matching GO IDs
for (name in names(comb_sub_list)) {
  matplot(x_list[[name]], comb_sub_list[[name]], type="l", col=cols, lwd=2, lty=1,
          xlab="n gene pairs", ylab="Cumulative % matching GO IDs",
          main=name, ylim=range(comb_list[[1]]))
}

# Bottom row: correlations
for (name in names(corr_sub_list)) {
  matplot(x_list[[name]], corr_sub_list[[name]], type="l", col=cols, lwd=2, lty=1,
          xlab="n gene pairs", ylab="Correlation", ylim=c(0,1),
          axes = FALSE)
  axis(1)
  axis(2, at=seq(0,1, by = 0.1))
  box()
  
  abline(h = c(0.7, 0.8, 0.9), col = "grey60", lty = 2, lwd = 1.5)
}

# ===================== Add shared legend =====================
par(xpd=NA)  # allow drawing outside the plotting area
legend("bottom", inset=-0.05, legend=methods, col=cols, lwd=2, lty=1,
       ncol=length(methods), cex=0.7)
