#' campp2_brca_1 - example data
#'
#' @description A campp2_brca_1 contains 10000 gene counts for 8 samples: 4 tumor (1
#' subtype) x 4 normal samples). Data are derived
#' from TCGA dataset.
#' @docType data
#' @usage data(campp2_brca_1)
#' @aliases campp2_brca_1
#' @return A data frame with 10000 rows (genes) and 8 columns (samples)
#' @format A data frame with 10000 rows (genes) and 8 columns (samples)
#'
"campp2_brca_1"


#' campp2_brca_2 - example data
#'
#' @description A campp2_brca_2 contains 10000 gene counts for 8 samples: 4 tumor (1
#' subtype) x 4 normal samples). Data are derived
#' from TCGA dataset.
#' @docType data
#' @usage data(campp2_brca_2)
#' @aliases campp2_brca_2
#' @return A data frame with 10000 rows (genes) and 8 columns (samples)
#' @format A data frame with 10000 rows (genes) and 8 columns (samples)
#'
"campp2_brca_2"


#' campp2_brca_1_meta - example data
#'
#' @description A metadata for campp2_brca_1 contains 9 variables for 8 samples in
#' campp2_brca_1 dataset. Data are derived from TCGA dataset.
#' @docType data
#' @usage data(campp2_brca_1_meta)
#' @aliases campp2_brca_1_meta
#' @return A data frame with 8 rows (samples) and 10 columns.
#' @format A data frame with 8 rows (samples) and 10 columns.
#'
"campp2_brca_1_meta"


#' campp2_brca_2_meta - example data
#'
#' @description A metadata for campp2_brca_2 contains 9 variables for 8 samples in
#' campp2_brca_2 dataset. Data are derived from TCGA dataset.
#' @docType data
#' @usage data(campp2_brca_2_meta)
#' @aliases campp2_brca_2_meta
#' @return A data frame with 8 rows (samples) and 10 columns.
#' @format A data frame with 8 rows (samples) and 10 columns.
#'
"campp2_brca_2_meta"


#' campp2_brca_1_NAs - example data
#'
#' @description A dataset with randomly introduced 1000 NA values into every sample's genes
#' (10000). This dataset is based on "campp2_brca_1". For more details, see the
#' vignette.
#' @docType data
#' @usage data(campp2_brca_1_NAs)
#' @aliases campp2_brca_1_NAs
#' @return A matrix with 10000 rows (genes) and 8 columns (samples).
#' @format A matrix with 10000 rows (genes) and 8 columns (samples).
#'
"campp2_brca_1_NAs"


#' campp2_brca_1_replacedNAs - example data
#'
#' @description A dataset with replaced NA values. This dataset is based on
#' "campp2_brca_1_NAs" dataset. For more details, see the vignette.
#' @docType data
#' @usage data(campp2_brca_1_replacedNAs)
#' @aliases campp2_brca_1_replacedNAs
#' @return A data frame with 10000 rows (genes) and 8 columns (samples).
#' @format A data frame with 10000 rows (genes) and 8 columns (samples).
#'
"campp2_brca_1_replacedNAs"


#' campp2_brca_1_zeroFix - example data
#'
#' @description A dataset with fixed zero values generated using "FixZeros" function. This
#' dataset is based on "campp2_brca_1_replacedNAs" data. For more details, see the vignette.
#' @docType data
#' @usage data(campp2_brca_1_zeroFix)
#' @aliases campp2_brca_1_zeroFix
#' @return A data frame with 9720 rows (genes) and 8 columns (samples).
#' @format A data frame with 9720 rows (genes) and 8 columns (samples).
#'
"campp2_brca_1_zeroFix"


#' campp2_brca_1_normalized - example data
#' @docType data
#' @usage data(campp2_brca_1_normalized)
#' @aliases campp2_brca_1_normalized
#' @description A dataset with normalized and transformed gene counts data generated using
#' "NormalizeData" function. This dataset is based on "campp2_brca_1_zeroFix" data.
#' For more details, see the vignette.
#' @return A Elist with 9626 rows (genes) and 8 columns (samples).
#' @format A Elist with 9626 rows (genes) and 8 columns (samples).
#'
"campp2_brca_1_normalized"


#' campp2_brca_1_batchCorrected - example data
#'
#' @description A dataset with batch corrected data generated using "BatchCorrect" function.
#' This dataset is based on "campp2_brca_1_normalized". For more details, see the
#' vignette.
#' @docType data
#' @usage data(campp2_brca_1_batchCorrected)
#' @aliases campp2_brca_1_batchCorrected
#' @return A matrix with 8919 rows (genes) and 8 columns (samples).
#' @format A matrix with 8919 rows (genes) and 8 columns (samples).
#'
"campp2_brca_1_batchCorrected"
#'
#'
#' campp2_brca_1_distributionsFit - example data
#' @docType data
#' @usage data(campp2_brca_1_distributionsFit)
#' @aliases campp2_brca_1_distributionsFit
#' @description A results from FitDistributions function providing a list of distribution
#' descriptions for first 10 features.
#' @return A list of 10 lists describing fitting the normal distribution.
#' @format A list of 10 lists describing fitting the normal distribution.
#'
"campp2_brca_1_distributionsFit"
#'

#' campp2_brca_1_batchCorrected_mean - example data
#' @docType data
#' @usage data(campp2_brca_1_batchCorrected_mean)
#' @aliases campp2_brca_1_batchCorrected_mean
#' @description This object includes information about an average counts and SD for
#' each feature across all the samples and for each sample group. Results are
#' provided in a form of a list of data frames.
#' Data were generated by running:
#' campp2_brca_1_batchCorrected_mean<-
#' meanCounts(campp2_brca_1_batchCorrected, campp2_brca_1_meta$diagnosis)
#' @return a list of 2 data frames (both 9048 x 7)
#' @format a list of 2 data frames (both 9048 x 7)
#'
"campp2_brca_1_batchCorrected_mean"
#'
#'

#' campp2_brca_1_DEA - example data
#' @docType data
#' @usage data(campp2_brca_1_DEA)
#' @aliases campp2_brca_1_DEA
#' @description Results from differential gene expression analysis provided by
#' RunDEA function using the provided example script:
#' campp2_brca_1_DEA<-RunDEA(data=campp2_brca_1, metadata=campp2_brca_1_meta,
#' group=campp2_brca_1_meta$subtype, prefix="test",
#' block=NULL, batch=campp2_brca_1_meta$age,
#' covarDEA = c("tumor_stage"), cutoff.logFC=1, cutoff.FDR=0.05)
#' @return a list of:
#' 1) a data frame (43 x 8) - re-formatted matrix of differential
#' expression/abundance results from limma.
#' 2) a list of original results (a list) from limma
#' 3) a character vector with unique feature names
#' 4) a design matrix (8x4)
#' 5) a contrast matrix (4x1)
#' @format a list of:
#' 1) a data frame (43 x 8) - re-formatted matrix of differential
#' expression/abundance results from limma.
#' 2) a list of original results (a list) from limma
#' 3) a character vector with unique feature names
#' 4) a design matrix (8x4)
#' 5) a contrast matrix (4x1)
#'
"campp2_brca_1_DEA"


#' campp2_brca_1_DEA_HUGO - example data
#' @docType data
#' @usage data(campp2_brca_1_DEA_HUGO)
#' @aliases campp2_brca_1_DEA_HUGO
#' @description results from differential gene expression analysis provided by
#' RunDEA, filtering on group comparisons and applying AddGeneName:
#' step1:
#' campp2_brca_1_DEA<-RunDEA(data=campp2_brca_1, metadata=campp2_brca_1_meta,
#' group=campp2_brca_1_meta$subtype, prefix="test",
#' block=NULL, batch=campp2_brca_1_meta$age,
#' covarDEA = c("tumor_stage"), cutoff.logFC=1, cutoff.FDR=0.05)
#' step2:
#' campp2_brca_1_DEA_HUGO<-subset(campp2_brca_1_DEA, grepl("healthy",
#' comparison, fixed = TRUE))
#' step3:
#' campp2_brca_1_DEA_HUGO<-AddGeneName(campp2_brca_1_DEA_HUGO,ensembl.version)
#' @return a data frame (43 x 9)
#' @format a data frame (43 x 9)
#'
"campp2_brca_1_DEA_HUGO"
#'
#'
#' campp2_brca_1_DEA_out - example data
#'
#' @description results from differential gene expression analysis provided by
#' ExportDEA function using the provided example script:
#' campp2_brca_1_DEA_out<-ExportDEA(res.DEA = campp2_brca_1_DEA$res.DEA,
#'                                  prefix="test")
#' @docType data
#' @usage data(campp2_brca_1_DEA_out)
#' @aliases campp2_brca_1_DEA_out
#' @return a list of:
#' 1) a data frame (43 x 8) - re-formatted matrix of differential
#' expression/abundance results from limma.
#' 2) a list of original results (a list) from limma
#' 3) a character vector with unique feature names
#' 4) a design matrix
#' 5) a contrast matrix
#' @format a list of:
#' 1) a data frame (43 x 8) - re-formatted matrix of differential
#' expression/abundance results from limma.
#' 2) a list of original results (a list) from limma
#' 3) a character vector with unique feature names
#' 4) a design matrix
#' 5) a contrast matrix
#'
"campp2_brca_1_DEA_out"
#'
#'
#' campp2_brca_1_LASSO - example data
#'
#' @description This object represents an output from runLASSO function. The list
#' includes a dataframe of feature names detected by both DEA and EN/LASSO/Ridge
#' regression. Another 2 sub-objects from the list for cross validation error
#' and roc. are empty because the input data
#' is not large enough for the computation of those metrics.
#' Data were generated by running because LASSO can be run on data where there
#' is at least 8 samples per group. Here we create a large dataset:
#' campp2_test_data_LASSO<-cbind(campp2_brca_1,campp2_brca_2)
#' campp2_test_metadata_LASSO<-rbind(campp2_brca_1_meta, campp2_brca_2_meta)
#' campp2_test_data_LASSO_replaceNAs<-ReplaceNAs(data=campp2_test_data_LASSO)
#' campp2_test_data_LASSO_zeroFix<-FixZeros(data=campp2_test_data_LASSO_replaceNAs,group=campp2_test_metadata_LASSO$diagnosis, remove.sparse.features=TRUE)
#' campp2_test_data_LASSO_normalized<-NormalizeData(data=campp2_test_data_LASSO_zeroFix,group=campp2_test_metadata_LASSO$diagnosis,standardize="TMM",transform="voom",technology="seq")
#' campp2_test_data_LASSO_batchCorrected<-BatchCorrect(data=campp2_test_data_LASSO_normalized,batch=campp2_test_metadata_LASSO$tumor_stage,group=campp2_test_metadata_LASSO$diagnosis,technology="seq")
#' run lasso on a large dataset
#' campp2_brca_1_LASSO<-runLASSO(data=campp2_test_data_LASSO_batchCorrected, group=campp2_test_metadata_LASSO$diagnosis, alpha=0.5, min.coef=0, nfolds=10, prefix="test")
#' @docType data
#' @usage data(campp2_brca_1_LASSO)
#' @aliases campp2_brca_1_LASSO
#' @return #' a list of 3 data frame:
#' 1) (42x1) and 2 empty variables
#' 2) NA
#' 3) NA
#' @format
#' a list of 3 data frame:
#' 1) (42x1) and 2 empty variables
#' 2) NA
#' 3) NA
#'
"campp2_brca_1_LASSO"


#' campp2_brca_1_DEA_LASSO_consensus - example data
#'
#' @description This object represents a data frame output from
#' campp2_brca_1_DEA_LASSO_consensus function. The table includes intersect
#' of the DEA and EN/LASSO/Ridge regression results. The format of the data
#' corresponds to the DEA (limma) output.
#' Data were generated by running:
#' ##here we prepare DEA results compatible (group="diagnosis) with
#' EN/LASS/Ridge regression results:
#' campp2_brca_1_DEA_diagnosis<-RunDEA(data=campp2_brca_1_normalized,
#' metadata=campp2_brca_1_meta,
#' group=campp2_brca_1_meta$diagnosis, prefix="test", batch=campp2_brca_1_meta$age,
#' covarDEA = c("tumor_stage"), cutoff.logFC=1, cutoff.FDR=0.05)
#' ##here we run RunDEA_LASSO_consensus function:
#' campp2_brca_1_DEA_LASSO_consensus <- RunDEA_LASSO_consensus(
#' campp2_brca_1_DEA_diagnosis$DEA.out,
#' campp2_brca_1_LASSO, group=campp2_brca_1_meta$diagnosis,
#' viridis.palette="turbo", "test")
#' @docType data
#' @usage data(campp2_brca_1_DEA_LASSO_consensus)
#' @aliases campp2_brca_1_DEA_LASSO_consensus
#' @return a data frame (7x8)
#' @format a data frame (7x8)
"campp2_brca_1_DEA_LASSO_consensus"



#' campp2_brca_1_forest_features.rda - example data
#'
#' @description This object represents output from the
#' ForestFeatures function which is used to fit random forest model
#' to the data and to perform feature selection using random forest.
#' Data were generated by running:
#' campp2_brca_1_forest_features <- ForestFeatures(seed = 123,
#' data = campp2_brca_1_batchCorrected,
#' group = campp2_brca_1_meta$diagnosis,
#' validation = TRUE,
#' test.train.ratio = 0.25,
#' num.trees.init = 5000,
#' num.trees.iterat = 2000)
#' @docType data
#' @usage data(campp2_brca_1_forest_features)
#' @aliases campp2_brca_1_forest_features
#' @return a list containing four elements
#' @format a list containing four elements
#'
"campp2_brca_1_forest_features"


#' campp2_brca_1_rf_apply.rda - example data
#'
#' @description This object represents output from the
#' RFApply function which is used to run the ForestFeatures function
#' across seeds.
#' Data were generated by running:
#' campp2_brca_1_rf_apply <-
#' RFApply(data = campp2_brca_1_batchCorrected,
#' group = campp2_brca_1_meta$diagnosis,
#' validation = TRUE,
#' test.train.ratio = 0.25,
#' num.trees.init = 5000,
#' num.trees.iterat = 2000)
#' @docType data
#' @usage data(campp2_brca_1_rf_apply)
#' @aliases campp2_brca_1_rf_apply
#' @format a list containing eight elements
#' @return a list containing eight elements
#'
"campp2_brca_1_rf_apply"


#' campp2_brca_1_run_rf.rda - example data
#'
#' @description This object represents output from the
#' RunRF function which is used to run the RFApply function.
#' Data were generated by running:
#' campp2_brca_1_run_rf <-
#' RunRF(data = campp2_brca_1_batchCorrected,
#' group = campp2_brca_1_meta$diagnosis,
#' split.size = 5,
#' test.train.ratio = 0.25,
#' num.trees.init = 5000,
#' num.trees.iterat = 2000)
#' @docType data
#' @usage data(campp2_brca_1_run_rf)
#' @aliases campp2_brca_1_run_rf
#' @return a list containing two elements
#' @format a list containing two elements
#'
"campp2_brca_1_run_rf"
