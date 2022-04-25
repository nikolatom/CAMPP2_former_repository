#' @title Apply DEAFeature
#' @description Apply differential expression/abundance analysis to all comparisons at once.
#' @param contrast an array containing contrasts between groups of interest.
#' @param data Gene count matrix (gene IDs as row names and sample IDs as columns). It's recommended to import gene counts using function "import_counts".
#' @param design a design matrix with all comparisons.
#' @param coLFC a number specifying the cutoff for LogFC.
#' @param coFDR a number specifying the cutoff for FDR.
#' @param block a vector of patient IDs if blocking is provided.
#' @param vector a vector of patient IDs; TRUE/FALSE statement specifying output format. If TRUE the function return a vector of feature IDs only.
#' @export
#' @import sva
#' @import limma
#' @seealso
#' @return DA/DE features
#' @examples \dontrun{
#' ...
#' }

DEAFeatureApply <- function(contrast, data, design, coLFC, coFDR, block=NULL, vector=FALSE) {
    features.l <- apply(contrast, 2, function(x) DEAFeature(x, data, design, coLFC, coFDR, block))

    if(vector == TRUE) {
        features <- do.call(rbind, lapply(features.l, function(x) do.call(rbind, x)))
        features <- unique(do.call(rbind, strsplit(rownames(features), "[.]"))[,2])
        return(features)
    }
    else {
        return(features.l)
    }
}
