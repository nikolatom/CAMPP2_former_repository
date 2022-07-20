#' @title Run differential expression/abundance analysis
#' @descriptionA A function for running differential expression/abundance analysis on a matrix data sample using limma.
#' @param data A raw gene count matrix from seq, array, ms or other technology (with gene IDs as row names and sample IDs as columns). It's recommended to import gene counts using function "import_counts".
#' @param technology a string vector of length 1 defining technology used for generating the data. Allowed types are: 'array', 'seq', 'ms' or 'other'.
#' @param batch A list of batches in the input samples
#' @param covarDEA Covariates to include in the analysis If multiple of these, they should be specified as a character vector.
#' @param group A string vector of specific sample groups derived from metadata column (e.g. diagnosis)
#' @param cutoff.logFC A cutoff value for the logarithmic Fold Change for each data sample (ratio of changes in expression data)
#' @param cutoff.FDR The false discovery rate for each data sample (the corrected p-value)
#' @param prefix a character defining a prefix of output file names.
#' @param block A vector or factor specifying a blocking variable. The block must be of same length as data and contain 2 or more options. For 2 datasets, the block can be defined as a vector of the two seperate blocks.
#' @export
#' @import limma
#' @import sva
#' @seealso
#' @return a matrix of differential expression/abundance results (containing information on gene regulation and significance of comparison)
#' @examples \dontrun{
#' ...
#' }


RunDEA <- function(data, technology, batch, covarDEA, group, cutoff.logFC, cutoff.FDR, prefix, block) {

    if (!(technology) %in% c("seq", "array", "ms", "other")) {
        stop("Defined technology is not supported.")
    }

    # Make design matrix
    if(is.null(covarDEA)) {
        if (is.null(batch)) {
            design <- model.matrix(~0+group)
            out.name <- "_DE"
        } else if (length(batch) != ncol(data)) {
            stop("Batch correction selected, but batches column does not match the samples!")
        } else {
            design <- model.matrix(~0+group+batch)
            out.name <- "_databatch_DE"
        }
    } else {
        if (is.null(batch)) {
            design.str <- "model.matrix(~0+group"
            out.name <- "_DE"
        }

        s <- lapply(split(as.matrix(df), col(df)), factor)
        my.names <- paste0("", colnames(df))
        my.names <- paste0(my.names, collapse = "+")
        design <- eval(parse(text=paste0(design.str,"+",my.names,")")))
    }

    # Making group contrasts
    combinations <- data.frame(t(combn(paste0("group",levels(as.factor(group))), 2)))
    combinations$contr <- apply(combinations[,colnames(combinations)], 1, paste, collapse = "-")
    contrast.matrix <- makeContrasts(contrasts=combinations$contr,levels=as.character(colnames(design)))

    # Apply DEA to all comparisons
    res.DEA <- DEAFeatureApply(contrast.matrix, data, design, cutoff.logFC, cutoff.FDR, block)

    # Write results out as .txt file
    if (!is.null(res.DEA)) {
        DEA.out <- ExportDEA(res.DEA, paste0(prefix, out.name))
        rownames(DEA.out) <- NULL
        res.DEA.names <- unique(DEA.out$name)
    } else {
        cat("No signficant DEA hits found. Check your cut-off for differential expression analysis, it may be that these are too stringent.")
    }

    if (technology[1] == "seq") {
        cnames <- colnames(data$E)
        data <- data.frame(data$E)
        colnames(data) <- cnames
    }

    return(list("DEA.out"=DEA.out,"res.DEA"=res.DEA,"res.DEA.names"=res.DEA.names))
}
