#' @title Differential expression/abundance features
#' @description Function for returning DEA features for DEAFeatureApply function.
#' @param contrast an array containing contrasts between groups of interest.
#' @param data Gene count matrix (gene IDs as row names and sample IDs as columns). It's recommended to import gene counts using function "import_counts".
#' @param design a design matrix with all comparisons.
#' @param coLFC a number specifying the cutoff for LogFC.
#' @param coFDR a number specifying the cutoff for FDR.
#' @param block a vector of patient IDs if blocking is provided.
#' @param vector a vector of patient IDs; TRUE/FALSE statement specifying output format, if TRUE the function return a vector of feature IDs only
#' @export
#' @import sva
#' @import limma
#' @seealso
#' @return DA/DE features for DAFeatureApply function
#' @examples \dontrun{
#' ...
#' }


DEAFeature <- function(contrast, data, design, coLFC, coFDR, block=NULL) {
    if(is.null(block)) {
        fit3 <- eBayes(contrasts.fit(lmFit(data, design), contrast))
    }
    else {
        corfit <- duplicateCorrelation(data, design, block=block)
        fit3 <- eBayes(contrasts.fit(lmFit(data, design, block = block, correlation=corfit$consensus), contrast))
    }
    tt <- topTable(fit3, coef=1, adjust='fdr', number=nrow(data))

    up <- tt[tt$logFC >= coLFC & tt$adj.P.Val < coFDR, ]
    down <- tt[tt$logFC <= -coLFC & tt$adj.P.Val < coFDR, ]

    up$name <- rownames(up)
    down$name <- rownames(down)

    up$dir <- rep("up", nrow(up))
    down$dir <- rep("down", nrow(down))

    final <- list(up, down)
    return(final)
}
