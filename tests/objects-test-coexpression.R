# Load packages
source("R/packages.R")
# Load functions
source("R/coexpression-functions.R")

# Make smaller toy data
raw_counts <- read.csv("data/k12_counts.csv", header=TRUE, row.names = 1) 
test_exprmat <- filterLowExpression(raw_counts[1:100,])

counts_sorted <- sort_expression(t(test_exprmat))
z_counts <- cmultRepl(counts_sorted,method = "CZM",output="p-counts") 
comp <- get_compositions(z_counts)
clr <- get_pcc(comp)

pr <- get_prop(counts_sorted)
pr_z <- get_prop(z_counts)

test_mats <- list()
test_mats$CLR = clr
test_mats$propr = pr
test_mats$propr_z = pr_z

saveRDS(test_mats,"tests/test_mats.rds")

test_gp <- list()
test_gp$CLR = get_genepairs(clr)
test_gp$propr = get_genepairs(pr)
test_gp$propr_z = get_genepairs(pr_z)

saveRDS(test_gp,"tests/test_gp.rds")
