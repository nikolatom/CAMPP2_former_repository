#' @title Apply DEA to all group comparisons
#' @description a function for applying differential expression/abundance
#' analysis to all group comparisons at once. The analysis is usually based on
#' normalized feature counts (e.g., genes), design matrix and contrast matrix.
#' Cutoffs for logFC and FDR; and blocking variable are optional. Using
#' parameter vector=TRUE, only features' IDs are output.
#' @param data a matrix of (transformed and normalized) feature counts from
#' "seq", "array", "ms" or "other" technology (with feature IDs as row names
#' and sample IDs as columns).
#' @param design.matrix a design matrix based on samples metadata (groups, covariates,
#' etc.).
#' @param contrast.matrix an array containing contrast.matrix between groups of
#' interest.
#' @param cutoff.logFC A cutoff value for the logarithmic fold change applied
#' to each feature. Default = 1.
#' @param cutoff.FDR a cutoff value for the false discovery rate (corrected
#' p-value) applied to each feature. Default = 0.01.
#' @param block a factor specifying blocking variable. The block must be of
#' the same length as data and contain 2 or more options (e.g. could be
#' represented by a column from a metadata file). Default = NULL.
#' @param vector a boolean value. If TRUE, only a vector of unique features IDs
#' will be output. Default = FALSE.
#' @export
#' @import sva
#' @import limma
#' @return a list (a vector in case vector=TRUE) of DEA features for all
#' comparisons
#' @examples {
#' DEA_all_comparisons <- DEAFeatureApply(data = campp2_brca_1_normalized,
#' design.matrix = campp2_brca_1_DEA$DEA.design.matrix, contrast.matrix =
#' campp2_brca_1_DEA$DEA.contrast.matrix, cutoff.logFC =1, cutoff.FDR =0.05,
#' block = NULL, vector = FALSE)
#' }

DEAFeatureApply <- function(data, design.matrix, contrast.matrix, cutoff.logFC, cutoff.FDR, block, vector=FALSE) {
    DEA.features.list <- apply(contrast.matrix, 2, function(x) DEAFeature(x, data, design.matrix, cutoff.logFC, cutoff.FDR, block))
    if(vector == TRUE) {
        DEA.features.vector <- do.call(rbind, lapply(DEA.features.list, function(x) do.call(rbind, x)))
        DEA.features.vector <- unique(do.call(rbind, strsplit(rownames(features), "[.]"))[,2])
        return(DEA.features.vector)
    } else {
        return(DEA.features.list)
    }
}
