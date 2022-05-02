#' @title Differentially abundant features - DA analysis
#' @description Function for returning DA/DE features for DAFeatureApply function
#' @param my.contrast a contrast between groups of interest
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
#' @return DA/DE features for DAFeatureApply function
#' @examples \dontrun{
#' ...
#' }


DAFeature <- function(my.contrast, my.data, my.design, coLFC, coFDR, my.block=NULL) {
    if(is.null(my.block)) {
        fit3 <- eBayes(contrasts.fit(lmFit(my.data, my.design), my.contrast))
    }
    else {
        corfit <- duplicateCorrelation(my.data, my.design, block=my.block)
        fit3 <- eBayes(contrasts.fit(lmFit(my.data, my.design, block = my.block, correlation=corfit$consensus), my.contrast))
    }
    tt <- topTable(fit3, coef=1, adjust='fdr', number=nrow(my.data))

    up <- tt[tt$logFC >= coLFC & tt$adj.P.Val < coFDR, ]
    down <- tt[tt$logFC <= -coLFC & tt$adj.P.Val < coFDR, ]

    up$name <- rownames(up)
    down$name <- rownames(down)

    up$dir <- rep("up", nrow(up))
    down$dir <- rep("down", nrow(down))

    final <- list(up, down)
    return(final)
}
