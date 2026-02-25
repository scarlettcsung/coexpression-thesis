#!/usr/bin/env Rscript

# Usage Rscript tests/test-coexpression.R

# Load packages
suppressMessages(source("R/packages.R"))
# Load functions
suppressMessages(source("R/coexpression-functions.R"))

librarian::shelf(crayon)

# ================= OBJECTS FOR TESTS =================

# Load test toy ata
raw_counts <- read.csv("data/k12_counts.csv", header=TRUE, row.names = 1) 
test_exprmat <- filterLowExpression(raw_counts[1:100,])
test_mats <- readRDS("tests/test_mats.rds")
test_gp <- readRDS("tests/test_gp.rds")

# Test gene pair function objects

# Toy correlation matrix
corr_mat <- matrix(
  c(
    1.0, 0.1, 0.2,  # Gene1
    0.1, 1.0, 0.3,  # Gene2
    0.2, 0.3, 1.0   # Gene3
  ),
  nrow = 3,
  ncol = 3,
  byrow = TRUE
)

# Add gene names
rownames(corr_mat) <- paste0("Gene", 1:3)
colnames(corr_mat) <- paste0("Gene", 1:3)

toy_genepairs <- data.frame(
  gene1 = c("Gene3", "Gene3", "Gene2"),
  gene2 = c("Gene2", "Gene1", "Gene1"),
  correlation = c(0.3, 0.2, 0.1),
  stringsAsFactors = FALSE
)

toy_genepairs <- as.data.frame(toy_genepairs)

# ================= TESTS =================
# NOTE: test_that is not used here as scripts are not in a package

# test get_genepairs()
gp <- get_genepairs(corr_mat)
# Remove row names, as they are irrelevant to result
rownames(gp) <- NULL
rownames(toy_genepairs) <- NULL

message("Testing get_genepairs")
if (!identical(gp, toy_genepairs)) {
  cat(red("FAIL: get_genepairs \n"))
} else {
  cat(green("PASS: get_genepairs \n"))
}

# test coexpr_analysis(). 
# COCOLOCUS methods were already validated against density plots, 
# so we are only looking at CLR and propr methods

co_expr <- suppressMessages(coexpr_analysis(test_exprmat,2))
ce_mat_clr <- co_expr$matrix$CLR
ce_gp_clr <- co_expr$genepairs$CLR

message("Testing coexpr_analysis")
if (!identical(ce_mat_clr, test_mats$CLR)) {
  cat(red("FAIL: coexpr_analysis CLR matrix \n"))
} else {
  cat(green("PASS: coexpr_analysis CLR matrix \n"))
}

if (!identical(ce_gp_clr, test_gp$CLR)) {
  cat(red("FAIL: coexpr_analysis CLR gene pairs\n"))
} else {
  cat(green("PASS: coexpr_analysis CLR gene pairs\n"))
}

ce_mat_pr <- co_expr$matrix$propr
ce_gp_pr <- co_expr$genepairs$propr

if (!identical(ce_mat_pr, test_mats$propr)) {
  cat(red("FAIL: coexpr_analysis propr matrix\n"))
} else {
  cat(green("PASS: coexpr_analysis propr matrix\n"))
}

if (!identical(ce_gp_pr, test_gp$propr)) {
  cat(red("FAIL: coexpr_analysis propr gene pairs\n"))
} else {
  cat(green("PASS: coexpr_analysis propr gene pairs\n"))
}

ce_mat_prz <- co_expr$matrix$propr_z
ce_gp_prz <- co_expr$genepairs$propr_z

if (!identical(ce_mat_prz, test_mats$propr_z)) {
  cat(red("FAIL: coexpr_analysis propr matrix\n"))
} else {
  cat(green("PASS: coexpr_analysis propr matrix\n"))
}

if (!identical(ce_gp_prz, test_gp$propr_z)) {
  cat(red("FAIL: coexpr_analysis propr gene pairs\n"))
} else {
  cat(green("PASS: coexpr_analysis propr gene pairs\n"))
}
