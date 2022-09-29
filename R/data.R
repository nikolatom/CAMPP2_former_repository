#' campp2_brca_1 - example data
#'
#' A campp2_brca_1 contains 10000 gene counts for 30 samples: 20 tumor (4
#' subtypes, each subtype has 5 samples) x 10 normal samples). Data are derived
#' from TCGA dataset.
#'
#' @format A data frame with 10000 rows (genes) and 30 columns (samples):
#' \describe{
#' }
#'
"campp2_brca_1"


#' campp2_brca_2 - example data
#'
#' A campp2_brca_2 contains 10000 gene counts for 30 samples: 20 tumor (4
#' subtypes, each subtype has 5 samples) x 10 normal samples). Data are derived
#' from TCGA dataset.
#'
#' @format A data frame with 10000 rows (genes) and 30 columns (samples):
#' \describe{
#' }
#'
"campp2_brca_2"


#' campp2_brca_1_meta - example data
#'
#' A metadata for campp2_brca_1 contains 9 variables for 30 samples in
#' campp2_brca_1 dataset. Data are derived from TCGA dataset.
#'
#' @format A data frame with 30 rows (samples) and 10 columns.
#' \describe{
#' }
#'
"campp2_brca_1_meta"


#' campp2_brca_2_meta - example data
#'
#' A metadata for campp2_brca_2 contains 9 variables for 30 samples in
#' campp2_brca_2 dataset. Data are derived from TCGA dataset.
#'
#' @format A data frame with 30 rows (samples) and 10 columns.
#' \describe{
#' }
#'
"campp2_brca_2_meta"


#' campp2_brca_1_NAs - example data
#'
#' A dataset with randomly introduced 1000 NA values into every sample's genes
#' (10000). This dataset is based on "campp2_brca_1". For more details, see the
#' vignette.
#'
#' @format A matrix with 10000 rows (genes) and 30 columns (samples).
#' \describe{
#' }
#'
"campp2_brca_1_NAs"


#' campp2_brca_1_replacedNAs - example data
#'
#' A dataset with replaced NA values. This dataset is based on
#' "campp2_brca_1_NAs" dataset. For more details, see the vignette.
#'
#' @format A data frame with 10000 rows (genes) and 30 columns (samples).
#' \describe{
#' }
#'
"campp2_brca_1_replacedNAs"


#' campp2_brca_1_zeroFix - example data
#'
#' A dataset with fixed zero values generated using "FixZeros" function. This
#' dataset is based on "campp2_brca_1_replacedNAs" data. For more details, see the vignette.
#'
#' @format A data frame with 9720 rows (genes) and 30 columns (samples).
#' \describe{
#' }
#'
"campp2_brca_1_zeroFix"


#' campp2_brca_1_normalized - example data
#'
#' A dataset with normalized and transformed gene counts data generated using
#' "NormalizeData" function. This dataset is based on "campp2_brca_1_zeroFix" data.
#' For more details, see the vignette.
#'
#' @format A Elist with 9046 rows (genes) and 30 columns (samples).
#' \describe{
#' }
#'
"campp2_brca_1_normalized"


#' campp2_brca_1_batchCorrected - example data
#'
#' A dataset with batch corrected data generated using "BatchCorrect" function.
#' This dataset is based on "campp2_brca_1_normalized". For more details, see the
#' vignette.
#'
#' @format A matrix with 9046 rows (genes) and 30 columns (samples).
#' \describe{
#' }
#'
"campp2_brca_1_batchCorrected"


#' campp2_brca_1_DEA - example data
#'
#' A results from differential gene expression analysis provided by
#' RunDEA function using the provided example script:
#' campp2_brca_1_DEA<-RunDEA(data=campp2_brca_1, metadata=campp2_brca_1_meta,
#' group=campp2_brca_1_meta$subtype, prefix="test",
#' block=campp2_brca_1_meta$subtype, batch=campp2_brca_1_meta$age,
#' covarDEA = c("tumor_stage"), cutoff.logFC=1, cutoff.FDR=0.01)
#' @format a list of:
#' 1) a data frame (2670 x 8) - re-formatted matrix of differential
#' expression/abundance results from limma.
#' 2) a list of original results (a list) from limma
#' 3) a character vector with unique feature names
#' 4) a design matrix
#' 5) a contrast matrix
#' \describe{
#' }
#'
"campp2_brca_1_DEA"


#' campp2_brca_1_DEA_HUGO - example data
#'
#' A results from differential gene expression analysis provided by
#' RunDEA, filtering on group comparisons and applying AddGeneName:
#' step1:
#' campp2_brca_1_DEA<-RunDEA(data=campp2_brca_1, metadata=campp2_brca_1_meta,
#' group=campp2_brca_1_meta$subtype, prefix="test",
#' block=campp2_brca_1_meta$subtype, batch=campp2_brca_1_meta$age,
#' covarDEA = c("tumor_stage"), cutoff.logFC=1, cutoff.FDR=0.01)
#' step2:
#' campp2_brca_1_DEA_HUGO<-subset(campp2_brca_1_DEA, grepl("healthy",
#' comparison, fixed = TRUE))
#' step3:
#' campp2_brca_1_DEA_HUGO<-AddGeneName(campp2_brca_1_DEA_HUGO,ensembl.version)
#' @format
#' a data frame (2054 x 9)
#' \describe{
#' }
#'
"campp2_brca_1_DEA_HUGO"


#' campp2_brca_1_DEA_HUGO_features_per_group - example data
#'
#' This object represents a list of the features characteristic for each group
#' and was created based on campp2_brca_1_DEA_HUGO (described above) following
#' these steps:
#' control.group="healthy"
#' groups.full.list <- split.data.frame(campp2_brca_1_DEA_HUGO, campp2_brca_1_DEA_HUGO$comparison)
#' campp2_brca_1_DEA_HUGO_features_per_group=list()
#' for (features in groups.full.list){
#'     campp2_brca_1_DEA_HUGO_features_per_group <- append(campp2_brca_1_DEA_HUGO_features_per_group,list(features$name))
#' }
#' names(campp2_brca_1_DEA_HUGO_features_per_group) <- gsub("-","",gsub(control.group,"",names(groups.full.list)))
#' @format
#' a list of 4 character vectors (1341, 395, 55, 263)
#' \describe{
#' }
#'
"campp2_brca_1_DEA_HUGO_features_per_group"

