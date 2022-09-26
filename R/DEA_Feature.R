#' @title DEA_feature - calculation of differential expression/abundance
#' based on one comparison
#' @description a function for applying differential expression/abundance
#' analysis using voom transform and limma. The analysis is based on feature
#' counts (e.g., genes), design matrix and contrast matrix (using only 1
#' contrast). Cutoffs for logFC and FDR; and blocking variable are optional.
#' @param contrast.matrix an array containing contrast.matrix between groups of
#' interest.
#' @param data a raw feature count matrix from "seq", "array", "ms" or "other"
#' technology (with feature IDs as row names and sample IDs as columns). It's
#' recommended to import feature counts using function "import_counts".
#' @param design.matrix a design matrix based on samples metadata (groups, covariates,
#' etc.).
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
#' @seealso
#' @return DEA features
#' @examples \dontrun{
#' DEAFeature(contrast.matrix = campp2_brca_1_DEA$DEA.contrast.matrix[,1],
#' data = campp2_brca_1, design.matrix = campp2_brca_1_DEA$DEA.design.matrix,
#' cutoff.logFC =1, cutoff.FDR =0.01, block = campp2_brca_1_meta$subtype)
#' }


DEAFeature <- function(contrast.matrix, data, design.matrix, cutoff.logFC = 1, cutoff.FDR = 0.01, block=NULL) {
    if(is.null(block)) {
        fit3 <- treat(contrasts.fit(lmFit(data, design.matrix), contrast.matrix))
    }
    else {
        corfit <- duplicateCorrelation(data, design.matrix, block=block)
        fit3 <- treat(contrasts.fit(lmFit(data, design.matrix, block = block, correlation=corfit$consensus), contrast.matrix))
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
