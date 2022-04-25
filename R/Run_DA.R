#' @title Differential expression/abundance analysis
#' @description Running differential gene expression analysis on a matrix data sample using limma
#' @param data Gene count matrix (gene IDs as row names and sample IDs as columns). It's recommended to import gene counts using function "import_counts".
#' @param metadata Samples' metadata table should be imported using function "import_metadata". Metadata must include exactly the same samples as gene counts (data1) and samples must be sorted similarly.
#' @param databatch Defines whether experimental batches are present
#' @param batch List of batches in the sample
#' @param covarD Covariates to include in the analysis
#' @param group A vector of integers specifying the group
#' @param logFC The logarithmic Fold Change (the ratio of changes in expression data)
#' @param FDR The false discovery rate calculated as the corrected p-value
#' @export
#' @import limma
#' @import sva
#' @seealso
#' @return
#' @examples \dontrun{
#' ...
#' }


RunDA <- function(data, metadata, technology, databatch, batch, covarD, group, logFC, FDR, prefix) {

    # Make design matrix
    if (databatch == "FALSE") {
        design.str <- "model.matrix(~0+group"
        out.name <- "_DE"
    } else if (length(batch) > 0) {
        design.str <- "model.matrix(~0+group+batch"
        out.name <- "_databatch_DE"
    } else {
        stop("Batch correction selected but no batches column found!")
    }


    if(is.null(covarD)) {
        design <- eval(parse(text=paste0(design.str, ")")))

    } else {
        if (length(covarD) == 1) {
            df <- data.frame(metadata[,colnames(metadata) %in% covarD])
            colnames(df) <- covarD
        } else {
            df <- metadata[,colnames(metadata) %in% covarD]
        }

        s <- lapply(split(as.matrix(df), col(df)), factor)
        my.names <- paste0("", colnames(df))
        list2env(setNames(s, my.names), envir=.GlobalEnv)
        my.names <- paste0(my.names, collapse = "+")
        design <- eval(parse(text=paste0(design.str,"+",my.names,")")))
    }



    # Making group contrasts
    combinations <- data.frame(t(combn(paste0("group", levels(group)), 2)))
    combinations$contr <- apply(combinations[,colnames(combinations)], 1, paste, collapse = "-")
    contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list(as.character(combinations$contr)),levels=list(design))))


    # Apply differential abundance analysis to all comparisons
    res.DE <- DAFeatureApply(contrast.matrix, data, design, logFC, FDR, NULL, FALSE)


    # Write results out as excel file
    if (!is.null(res.DE)) {
        DE.out <- TextOutput(res.DE, paste0(prefix, out.name))
        rownames(DE.out) <- NULL
        res.DE.names <- unique(DE.out$name)
    } else {
        cat("No signficant DE/DA hits found. Check output file from differential expression analysis. Check your cut-off for differential expression analysis, it may be that these are too stringent.")
    }

    if (technology[1] == "seq") {
        cnames <- colnames(data$E)
        data <- data.frame(data$E)
        colnames(data) <- cnames
    }

    return(list("DE.out"=DE.out,"res.DE"=res.DE,"res.DE.names"=res.DE.names))
}
