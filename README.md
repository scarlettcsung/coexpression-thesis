# MSc Thesis: Proportionality as a promising method to enhancing functional inference by co-expression analysis
## Abstract
Co-expression analysis is widely used for gene function inference. Conventional correlation calculation methods in co-expression analysis, however, rely on normalization strategies that have been found to introduce spurious correlation effects. A proposed explanation suggests that expression data is compositional and should thus be analyzed with consideration of the relative relationships of its components, such as genes, rather than using conventional correlation measures. In this report, different workflows incorporating different combinations of correlation calculations, proportionality measures, and data transformation methods are compared. The correlation distributions produced by the different workflows were investigated to identify signs of compositional effects. Correlation-based co-expression analysis methods are systematically compared to proportionality-based methods on their ability to accurately recover established functional and structural relationships between genes. This report shows that proportionality-based methods more effectively recover functional and structural relationships between genes in co-expression analysis, and thus demonstrate promise in improving the robustness of co-expression analysis in bacterial species. The resulting potential improvements in gene function prediction robustness could further expedite the identification of novel gene candidates for biological, agricultural, and medical applications. 

## Thesis Report
https://edepot.wur.nl/711344

## Folder Descriptions
- R/ : co-expression analysis script, functions, packages
- data/ : downloaded data. Raw counts are not included but sourced in raw_counts_sources.txt
- dependencies/ : packages installed from project
- evaluation/ : evaluation and plot generation scripts, only to be run from IDE
- reference/ : reference files
- results/ : co-expression analysis results
- tests/ : unit testing

## Usage
The following files can be run via. command line. 
- R/01-coexpression-analysis.R

  `Rscript R/01-coexpression-analysis.R {input_filepath} {dataset_name}`

- test/test-coexpression.R

  `Rscript test/test-coexpression.R`
  
- test/test-GO-match.R

  `Rscript test/test-GO-match.R`
