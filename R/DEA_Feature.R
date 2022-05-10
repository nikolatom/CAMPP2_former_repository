#' @title DEA_feature
#' @description A function for returning DEA features for DEAFeatureApply function.
#' @param contrast an array containing contrasts between groups of interest.
#' @param data A raw gene count matrix from seq, array, ms or other technology (with gene IDs as row names and sample IDs as columns). It's recommended to import gene counts using function "import_counts".
#' @param design a design matrix with all comparisons.
#' @param coLFC a number specifying the cutoff for LogFC.
#' @param coFDR a number specifying the cutoff for FDR.
#' @param block a vector or factor specifying a blocking variable
#' @param vector a vector of patient IDs; TRUE/FALSE statement specifying output format, if TRUE the function return a vector of feature IDs only
#' @export
#' @import sva
#' @import limma
#' @seealso
#' @return DEA features
#' @examples \dontrun{
#' ...
#' }


DEAFeature <- function(contrast, data, design, coLFC, coFDR, block=NULL) {
    if(is.null(block)) {
        fit3 <- eBayes(contrasts.fit(lmFit(data, design), contrast))
    }
    else {
        corfit <- dup.reglicateCorrelation(data, design, block=block)
        fit3 <- eBayes(contrasts.fit(lmFit(data, design, block = block, correlation=corfit$consensus), contrast))
    }
    DEA.table <- topTable(fit3, coef=1, adjust='fdr', number=nrow(data))

    up.reg <- DEA.table[DEA.table$logFC >= coLFC & DEA.table$adj.P.Val < coFDR, ]
    down.reg <- DEA.table[DEA.table$logFC <= -coLFC & DEA.table$adj.P.Val < coFDR, ]

    up.reg$name <- rownames(up.reg)
    down.reg$name <- rownames(down.reg)

    up.reg$dir <- rep("up.reg", nrow(up.reg))
    down.reg$dir <- rep("down.reg", nrow(down.reg))

    all.reg <- list(up.reg, down.reg)
    return(all.reg)
}
