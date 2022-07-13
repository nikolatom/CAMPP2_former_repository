#' @title DEA_feature
#' @description A function for returning DEA features for DEAFeatureApply function.
#' @param contrast an array containing contrasts between groups of interest.
#' @param data A raw gene count matrix from seq, array, ms or other technology (with gene IDs as row names and sample IDs as columns). It's recommended to import gene counts using function "import_counts".
#' @param design a design matrix with all comparisons.
#' @param cutoff.logFC a number specifying the cutoff for LogFC.
#' @param cutoff.FDR a number specifying the cutoff for FDR.
#' @param block A vector or factor specifying a blocking variable. The block must be of same length as data and contain 2 or more options. For 2 datasets, the block can be defined as a vector of the two seperate blocks.
#' @export
#' @import sva
#' @import edgeR
#' @import statmod
#' @seealso
#' @return DEA features
#' @examples \dontrun{
#' ...
#' }


DEAFeature <- function(contrast, data, design, cutoff.logFC, cutoff.FDR, block) {
    if(is.null(block)) {
        fit3 <- treat(contrasts.fit(lmFit(data, design), contrast))
    }
    else {
        corfit <- duplicateCorrelation(data, design, block=block)
        fit3 <- treat(contrasts.fit(lmFit(data, design, block = block, correlation=corfit$consensus), contrast))
    }
    DEA.table <- topTreat(fit3, coef=1, adjust='fdr', number=nrow(data))

    up.reg <- DEA.table[DEA.table$logFC >= cutoff.logFC & DEA.table$adj.P.Val < cutoff.FDR, ]
    down.reg <- DEA.table[DEA.table$logFC <= -cutoff.logFC & DEA.table$adj.P.Val < cutoff.FDR, ]

    up.reg$name <- rownames(up.reg)
    down.reg$name <- rownames(down.reg)

    up.reg$dir <- rep("up.reg", nrow(up.reg))
    down.reg$dir <- rep("down.reg", nrow(down.reg))

    all.reg <- list(up.reg, down.reg)

    return(all.reg)
}
