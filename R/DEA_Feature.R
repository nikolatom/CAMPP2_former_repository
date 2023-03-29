#' @title DEA_feature - calculation of differential expression/abundance
#' based on one comparison
#' @description a function for applying differential expression/abundance
#' analysis using limma. The analysis is usually based on normalized feature
#' counts (e.g., genes), design matrix and contrast matrix (using only 1
#' contrast). Cutoffs for logFC and FDR; and blocking variable are optional.
#' @param contrast.matrix an array containing contrast.matrix between groups of
#' interest.
#' @param data a matrix of (transformed and normalized) feature counts from
#' "seq", "array", "ms" or "other" technology (with feature IDs as row names
#' and sample IDs as columns).
#' @param design.matrix a design matrix based on samples metadata (groups,
#' covariates, etc.).
#' @param cutoff.logFC A cutoff value for the logarithmic fold change applied
#' to each feature. Default = 1.
#' @param cutoff.FDR a cutoff value for the false discovery rate (corrected
#' p-value) applied to each feature. Default = 0.01.
#' @param block a factor specifying blocking variable. The block must be of
#' the same length as data and contain 2 or more options (e.g. could be
#' represented by a column from a metadata file). Default = NULL.
#' @export
#' @import sva
#' @import edgeR
#' @import statmod
#' @import limma
#' @return a list of up and donw-regulated features in limma format
#' @examples {
#' DEA_one_comparison <-
#' DEAFeature(contrast.matrix = campp2_brca_1_DEA$DEA.contrast.matrix[,1],
#' data = campp2_brca_1_normalized,
#' design.matrix = campp2_brca_1_DEA$DEA.design.matrix,
#' cutoff.logFC =1, cutoff.FDR =0.05, block = NULL)
#' }


DEAFeature <- function(contrast.matrix, data, design.matrix, cutoff.logFC = 1, cutoff.FDR = 0.01, block=NULL) {
    if(is.null(block)) {
        fit3 <- treat(contrasts.fit(lmFit(data, design.matrix), contrast.matrix), lfc = cutoff.logFC)
    }
    else {
        corfit <- duplicateCorrelation(data, design.matrix, block=block)
        fit3 <- treat(contrasts.fit(lmFit(data, design.matrix, block = block, correlation=corfit$consensus), contrast.matrix), lfc = cutoff.logFC)
    }
    DEA.table <- topTreat(fit3, coef=1, adjust='fdr', number=nrow(data))

    up.reg <- DEA.table[DEA.table$logFC >= cutoff.logFC & DEA.table$adj.P.Val < cutoff.FDR, ]
    down.reg <- DEA.table[DEA.table$logFC <= -(cutoff.logFC) & DEA.table$adj.P.Val < cutoff.FDR, ]

    up.reg$name <- rownames(up.reg)
    down.reg$name <- rownames(down.reg)

    up.reg$dir <- rep("up.reg", nrow(up.reg))
    down.reg$dir <- rep("down.reg", nrow(down.reg))

    DEA.features.list <- list(up.reg, down.reg)

    return(DEA.features.list)
}
