#!/usr/bin/env Rscript

# Usage Rscript tests/test-GO-match.R

suppressMessages(source("R/packages.R"))
suppressMessages(source("R/GO-overlap-functions.R"))

librarian::shelf(crayon)

# ================= CREATE TEST OBJECTS =================
go_ref <- data.frame(
  locus_tag = c("g1","g1","g1",
                "g2","g3","g4","g5","g6"),
  GO_ID = c("GO:0071239","GO:0006865","GO:0055085", #g1
            "GO:0071239", #g2
            "GO:0006865", #g3
            "GO:0071239", #g4
            "GO:0046373", #g5
            "GO:0046373"), #g6
  type = "process"
            
)

test_gp <- data.frame(
  gene1 =       c("g1","g1","g1","g2","g3","g4","g5","g6","g7"),
  gene2 =       c("g2","g3","g4","g5","g6","g7","g8","g1","g1"),
  correlation = c(0.9,0.89,0.87,0.7,0.6,0.6,0.5,0.4,0.3))

test_filtered_gp <- data.frame(
  gene1 =       c("g1","g1","g1","g2","g3","g6"),
  gene2 =       c("g2","g3","g4","g5","g6","g1"),
  correlation = c(0.9,0.89,0.87,0.7,0.6,0.4))

mapped_gp <- data.frame(
  gene1 =       c("g1","g1","g1","g2","g3","g6"),
  gene2 =       c("g2","g3","g4","g5","g6","g1"),
  GO_IDs_gene1 = c("GO:0071239;GO:0006865;GO:0055085",
                   "GO:0071239;GO:0006865;GO:0055085",
                   "GO:0071239;GO:0006865;GO:0055085",
                   "GO:0071239","GO:0006865","GO:0046373"),
  GO_IDs_gene2 = c("GO:0071239","GO:0006865","GO:0071239",
                   "GO:0046373","GO:0046373",
                   "GO:0071239;GO:0006865;GO:0055085"),
  correlation = c(0.9,0.89,0.87,0.7,0.6,0.4))


# ================= TESTS =================
# NOTE: test_that is not used here as scripts are not in a package

# test filter_process_go()------------------------------------------------------
filtered <- suppressMessages(filter_process_go(test_gp,go_ref))

# Remove row names, as they are irrelevant to result
rownames(filtered) <- NULL
rownames(test_filtered_gp) <- NULL

message("Testing filtered_process_go()")
if (!identical(filtered, test_filtered_gp)) {
  cat(red("FAIL: filter_process_go"))
} else {
  cat(green("PASS: filter_process_go"))
}

# test get_process_go_ids()-----------------------------------------------------
mapped <- suppressMessages(get_process_go_ids(filtered,go_ref,
                                              cor_values = filtered$correlation))

# Remove row names, as they are irrelevant to result
rownames(mapped) <- NULL
rownames(mapped_gp) <- NULL

message("Testing get_process_go_ids()")
if (!identical(mapped, mapped_gp)) {
  cat(red("FAIL: get_process_go_ids"))
} else {
  cat(green("PASS: get_process_go_ids"))
}
