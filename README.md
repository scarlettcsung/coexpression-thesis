# MSc Thesis: Proportionality as a promising method to enhancing functional inference by co-expression analysis
## Abstract
Co-expression analysis is widely used for gene function inference. Conventional correlation calculation methods in co-expression analysis, however, rely on normalization strategies that have been found to introduce spurious correlation effects. A proposed explanation suggests that expression data is compositional and should thus be analyzed with consideration of the relative relationships of its components, such as genes, rather than using conventional correlation measures. In this report, different workflows incorporating different combinations of correlation calculations, proportionality measures, and data transformation methods are compared. The correlation distributions produced by the different workflows were investigated to identify signs of compositional effects. Correlation-based co-expression analysis methods are systematically compared to proportionality-based methods on their ability to accurately recover established functional and structural relationships between genes. This report shows that proportionality-based methods more effectively recover functional and structural relationships between genes in co-expression analysis, and thus demonstrate promise in improving the robustness of co-expression analysis in bacterial species. The resulting potential improvements in gene function prediction robustness could further expedite the identification of novel gene candidates for biological, agricultural, and medical applications. 

## Project Structure
coexpression-thesis/
├── coexpression-thesis.Rproj 
├── data/bacillus_counts.csv 
├── data/BURK_counts.csv 
├── data/k12_counts.csv 
├── data/known_operon.download.txt 
├── data/raw_counts_sources 
├── dependencies/COBAIN.tar.gz 
├── evaluation/density_plots.R 
├── evaluation/eval_results/clustering/GO_terms_per_method.csv 
├── evaluation/eval_results/clustering/k12_clusters.rds 
├── evaluation/eval_results/clustering/k12_ego.rds 
├── evaluation/eval_results/clustering/k12_ego_simpl.rds 
├── evaluation/eval_results/clustering/k12_filtered_clusters.rds 
├── evaluation/eval_results/comp_benchmark_results.rds 
├── evaluation/eval_results/corr_distributions/binned_combined.rds 
├── evaluation/eval_results/corr_distributions/coexpression_medians.csv 
├── evaluation/eval_results/figures/density_coexpression.png 
├── evaluation/eval_results/permutation_testing/bac_background_means_10000perms.rds 
├── evaluation/eval_results/permutation_testing/bac_pvals_10000perms.rds 
├── evaluation/eval_results/permutation_testing/k12_background_means_10000perms.rds 
├── evaluation/eval_results/permutation_testing/k12_ko_background_means_10000perms.rds 
├── evaluation/eval_results/permutation_testing/k12_ko_pvals_10000perms.rds 
├── evaluation/eval_results/permutation_testing/k12_pvals_10000perms.rds 
├── evaluation/GO_enrichment.R 
├── evaluation/GO_matches.R 
├── evaluation/GOID_match_plots.R 
├── evaluation/KOID_match_plots.R 
├── evaluation/operon_matches.R 
├── evaluation/operon_plots.R 
├── evaluation/permutation_testing.R 
├── evaluation/runtime_memory.R 
├── evaluation/threshold_selection_plots.R 
├── R/01-coexpression-analysis.R 
├── R/coexpression-functions.R 
├── R/GO-overlap-functions.R 
├── R/KOID-overlap-functions.R 
├── R/packages.R 
├── R/permutation_testing_functions.R 
├── R/plotting-functions.R 
├── reference/derived/bacillus_genenames_sorted.rds 
├── reference/derived/burkholderia_genenames_sorted.rds 
├── reference/derived/go_bacillus.tsv 
├── reference/derived/go_burkholderia.tsv 
├── reference/derived/go_k12.tsv 
├── reference/derived/k12_gene_name_ref.tsv 
├── reference/derived/k12_genenames_sorted.rds 
├── reference/derived/k12_operon_reference.csv 
├── reference/raw/k12_known_operon.download.txt 
├── results/bacillus_filtered.rds 
├── results/burkholderia_filtered.rds 
├── results/k12_filtered.rds 
├── results_ignore/bacillus_genepairs.rds 
├── results_ignore/bacillus_GOmatches.rds 
├── results_ignore/bacillus_matrices.rds 
├── results_ignore/burkholderia_genepairs.rds 
├── results_ignore/burkholderia_GOmatches.rds 
├── results_ignore/burkholderia_matrices.rds 
├── results_ignore/k12_genepairs.rds 
├── results_ignore/k12_GOmatches.rds 
├── results_ignore/k12_komatches.rds 
├── results_ignore/k12_matrices.rds 
├── results_ignore/k12_shared.rds 
├── ScarletSung_MScThesis.pdf 
├── tests/objects-test-coexpression.R 
├── tests/test-coexpression.R 
├── tests/test-GO-match.R 
├── tests/test_gp.rds 
├── tests/test_mats.rds 

## Usage
The following files can be run via. command line. 
- R/01-coexpression-analysis.R

  `Rscript R/01-coexpression-analysis.R {input_filepath} {dataset_name}`

- test/test-coexpression.R

  `Rscript test/test-coexpression.R`
  
- test/test-GO-match.R

  `Rscript test/test-GO-match.R`
