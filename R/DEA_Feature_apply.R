#' @title Apply DEA_Feature
#' @description a function for applying differential expression/abundance
#' analysis to all group comparisons at once.
#' @param contrast.matrix an array containing contrast.matrix between groups of interest.
#' @param data a raw gene count matrix from seq, array, ms or other technology (with gene IDs as row names and sample IDs as columns). It's recommended to import gene counts using function "import_counts".
#' @param design a design matrix with all comparisons.
#' @param cutoff.logFC a number specifying the logFC cutoff for each feature.
#' @param cutoff.FDR a number specifying the cutoff for FDR.
#' @param block a factor (or a vector) specifying blocking variable. The block must be of
#' same length as data and contains 2 or more levels.
#' @param vector a boolean value. If TRUE, only a vector of features IDs will be output. Default = FALSE.
#' @export
#' @import sva
#' @import limma
#' @seealso
#' @return a list (alternatively a vector) of DEA features for all comparisons
#' @examples \dontrun{
#' ...
#' }

DEAFeatureApply <- function(data, design, contrast.matrix, cutoff.logFC, cutoff.FDR, block, vector=FALSE) {
    features.list <- apply(contrast.matrix, 2, function(x) DEAFeature(x, data, design, cutoff.logFC, cutoff.FDR, block))

    if(vector == TRUE) {
        features <- do.call(rbind, lapply(features.list, function(x) do.call(rbind, x)))
        features <- unique(do.call(rbind, strsplit(rownames(features), "[.]"))[,2])
        return(features)
    } else {
        return(features.list)
    }
}
