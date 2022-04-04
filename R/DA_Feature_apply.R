#' @title Apply DAFeature
#' @description Apply differential abundance analysis to all the comparisons at once
#' @param my.contrasts = a contrast between groups of interest
#' @param my.data a dataframe of expression/abundance counts
#' @param my.design a design matrix with all comparisons
#' @param my.coLFC a cutoff for logFC
#' @param my.coFDR a cutoff for FDR
#' @param my.block if blocking than a vector of patient IDs
#' @param my.vector a vector of patient IDs; TRUE/FALSE statement specifying output format, if TRUE the function return a vector of feature IDs only
#' @export
#' @import sva
#' @import limma
#' @seealso
#' @return DA/DE features
#' @examples \dontrun{
#' ...
#' }

DAFeatureApply <- function(my.contrasts, my.data, my.design, coLFC, coFDR, my.block=NULL, my.vector=FALSE) {
    my.features.l <- apply(my.contrasts, 2, function(x) DAFeature(x, my.data, my.design, coLFC, coFDR, my.block))
    print("test my.vector")
    print(my.vector)
    if(my.vector == TRUE) {
        my.features <- do.call(rbind, lapply(my.features.l, function(x) do.call(rbind, x)))
        my.features <- unique(do.call(rbind, strsplit(rownames(my.features), "[.]"))[,2])
        return(my.features)
    }
    else {
        return(my.features.l)
    }
}

