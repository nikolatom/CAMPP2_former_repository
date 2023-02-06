FROM bioconductor/bioconductor_docker:devel

WORKDIR /home/rstudio

RUN MAKEFLAGS="-j 8" R -e "library(BiocManager); BiocManager::install(version='devel', ask=FALSE); BiocManager::install(c('rio', 'factoextra', 'WGCNA', 'M3C', 'FactoMineR', 'statmod', 'edgeR', 'glmnet', 'sva', 'openxlsx', 'ggplot2', 'ggrepel', 'gridExtra', 'plyr', 'dplyr', 'data.table', 'viridis', 'squash', 'survcomp', 'scales', 'GO.db', 'fitdistrplus', 'impute', 'pcaMethods', 'pROC', 'VennDiagram', 'grid', 'ComplexHeatmap', 'SummarizedExperiment', 'mclust', 'multiMiR', 'biomaRt', 'doMC', 'nnet', 'zeallot', 'igraph', 'varSelRF', 'randomForest', 'caret', 'BiocStyle', 'knitr', 'rmarkdown'), ask=FALSE)"

