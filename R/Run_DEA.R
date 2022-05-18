#' @title Run differential expression/abundance analysis
#' @descriptionA A function for running differential expression/abundance analysis on a matrix data sample using limma.
#' @param data A raw gene count matrix from seq, array, ms or other technology (with gene IDs as row names and sample IDs as columns). It's recommended to import gene counts using function "import_counts".
#' @param metadata Samples' metadata table is recommended to be imported using function "import_metadata". Metadata must include exactly the same samples sorted in the same order as in a gene counts matrix (data).
#' @param technology a string vector of length 1 defining technology used for generating the data. Allowed types are: 'array', 'seq', 'ms' or 'other'.
#' @param databatch TRUE/FALSE value indicating activation/deactivation of batch correction. Boolean value is automatically set based on the definition of the "batches" parameter.
#' @param batch The batch covariate for each data sample, derived from metadata column.
#' @param covarDEA Covariates to include in the analysis If multiple of these, they should be specified as a character vector.
#' @param group An array of specific sample groups derived from metadata column (e.g. diagnosis)
#' @param logFC The logarithmic Fold Change for each data sample (ratio of changes in expression data)
#' @param FDR The false discovery rate for each data sample (the corrected p-value)
#' @param prefix a character defining the result folder name and prefix of output file names.
#' @export
#' @import limma
#' @import sva
#' @seealso
#' @return a matrix of differential expression/abundance results (containing information on gene regulation and significance of comparison)
#' @examples \dontrun{
#' ...
#' }


RunDEA <- function(data, metadata, technology, databatch, batch, covarDEA, group, logFC, FDR, prefix) {

    if (!(technology) %in% c("seq", "array", "ms", "other")) {
        stop("Defined technology is not supported.")
    }

    # Make design matrix
    if(is.null(covarDEA)) {
    if (databatch == "FALSE") {
        design <- model.matrix(~0+group)
        out.name <- "_DE"
    } else if (length(batch) != ncol(data)) {
        stop("Batch correction selected, but batches column does not match the samples!")
    } else {
        design <- model.matrix(~0+group+batch)
        out.name <- "_databatch_DE"
    }
    } else {
        if (databatch == "FALSE") {
            design.str <- "model.matrix(~0+group"
            out.name <- "_DE"
        }

        s <- lapply(split(as.matrix(df), col(df)), factor)
        my.names <- paste0("", colnames(df))
        my.names <- paste0(my.names, collapse = "+")
        design <- eval(parse(text=paste0(design.str,"+",my.names,")")))
    }


    # Making group contrasts
    combinations <- data.frame(t(combn(paste0("group", levels(group)), 2)))
    combinations$contr <- apply(combinations[,colnames(combinations)], 1, paste, collapse = "-")
    contrast.matrix <- makeContrasts(contrasts=combinations$contr,levels=as.character(colnames(design)))

    # Apply DEA to all comparisons
    res.DEA <- DEAFeatureApply(contrast.matrix, data, design, logFC, FDR, NULL, FALSE)


    # Write results out as .txt file
    if (!is.null(res.DEA)) {
        DEA.out <- TextOutput(res.DEA, paste0(prefix, out.name))
        rownames(DEA.out) <- NULL
        res.DEA.names <- unique(DEA.out$name)
    } else {
        cat("No signficant DE/DA hits found. Check your cut-off for differential expression analysis, it may be that these are too stringent.")
    }

    if (technology[1] == "seq") {
        cnames <- colnames(data$E)
        data <- data.frame(data$E)
        colnames(data) <- cnames
    }

    return(list("DEA.out"=DEA.out,"res.DEA"=res.DEA,"res.DEA.names"=res.DEA.names))
}
