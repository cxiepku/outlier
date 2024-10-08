# outlier
extreme outlier RNA expression

This repository includes the R code and test data for the manuscript "Patterns of extreme outlier RNA expression in population data reveal sporadic over-activation of genes with co-regulated modules in subsets of individuals". 

The R script "analyze_outlier.R" identifies outlier genes and correlated over outlier (OO) gene pairs. It takes a TPM matrix with gene annotation (example: "data/GTEx_brain_TPM.tsv"), a sample annotation file (example: "data/GTEx_brain_coldata.tsv"), and a list of genes with high sequence similarity paralogs to be excluded (example: "data/human_genes_wParalogs.list"; see the manuscript for further explanation) as input files. It also takes a few input parameters (see the script). It generates two output files. One includes the outlier genes (example: "data/GTEx_brain_outlier.tsv"; note that it also contains under outliers - UOs), and the other includes the correlated OO gene pairs (example: "data/GTEx_brain_outlier_pair.tsv"). It only requires two R packages: matrixStats and dplyr, and runs very fast. It has been tested under R 4.2.3, matrixStats 1.3.0, and dplyr 1.1.4, but the versions do not matter much.

The test data mentioned above is in the folder "data/", which is the GTEx brain data (see the manuscript for further information). The full data is in the supplemental materials of the manuscript, Edmond, or ENA (see the manuscript).
