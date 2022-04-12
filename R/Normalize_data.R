#' @title Normalization and Transformation
#' @description This function normalizes and transforms feature counts depending on data types and selected normalization/standardization methods.
#' @param data a dataframe with expression/abundance counts (genes as rows; samples as columns), N.B only a subset of variables should be input, not intended for the full expression matrix!
#' @param technology a string vector of length 1 (or two in case of 2 datasets) defining technology used for generating the data. Allowed types are: 'array', 'seq', 'ms' or 'other'.
#' @param group a dataframe with groups metadata for each sample (given from metadata table)
#' @param transform a string vector of length 1 (or two in case of 2 datasets) defining transformation type for each dataset ("log2", "logit", "log10"). voom transform is automatically applied on "seq" data.
#' @param standardize a string vector of length 1 (or two in case of 2 datasets) defining standardization method  ("mean" or "median"); By default, data from "seq" technology are normalized by "TMM"; "array" technology is standardized using "quantile"
#' @param data.original an original data1 (data.frame) including zero values
#' @export
#' @import fitdistrplus
#' @seealso
#' @return Elist of normalized, filtered and transformed (gene) counts data
#' @examples \dontrun{
#' ...
#' }
#'


NormalizeData <- function(data,data.original,group,transform,standardize,technology) {
    if (technology == "seq") {
        if (!is.null(data.original)) {
            data <- data.original
        }
        data <- DGEList(counts=data)
        design <- model.matrix(~0+group)
        keep <- filterByExpr(data, design)
        data <- data[keep,,keep.lib.sizes=FALSE]
        data <- calcNormFactors(data, method = "TMM")
        data <- voom(data, design, plot=TRUE)
        cat("\n Data will be filtered for lowly expressed variables, normalized and voom transformed.\n")


    } else if (technology %in% c("array", "ms", "other")) {
        if (transform == "log2") {
            data <- log2(data)
        } else if (transform == "logit") {
            data <- logit(data)
        } else if (transform == "log10"){
            data <- log10(data)
        } else {
            cat("\n transform is not specified for data, log transformation will NOT be performed.\n")
            data <- data
        }
        if (standardize == "mean") {
            data <- scale(data, scale = FALSE)
            NB <- " N.B This pipeline does not handle background correction of single-channel intensity data or within array normalization two-color intensity data. See limma manual section on array normalization for more on this. Data may be fully normalized with limma (R/Rstudio) or another software and the pipeline re-run."
            cat(paste0("\n technology = array and standardize = mean. Data will be mean centered.", NB, "\n"))

        } else if (standardize == "median") {
            rowmed <- apply(data,1,median)
            data <- data - rowmed
            NB <- " N.B This pipeline does not handle background correction of single-channel intensity data or within array normalization two-color intensity data. See limma manual section on array normalization for more on this. Data may be fully normalized with limma (R/Rstudio) or another software and the pipeline re-run."
            cat(paste0("\n technology = array and standardize = median. Data will be median centered.", NB, "\n"))
        } else if (!(standardize %in% c("mean", "median")) & technology == "array") {
            data <- normalizeBetweenArrays(data, method="quantile")
            NB <- " N.B This pipeline does not handle background correction of single-channel intensity data or within array normalization two-color intensity data. See limma manual section on array normalization for more on this. Data may be fully normalized with limma (R/Rstudio) or another software and the pipeline re-run."
            cat(paste0("\n technology = array. Data will be normalized on the quantiles.",NB, "\n"))
        } else {
            cat("\n- No standardization requested. If argument technology is 'array', data will be normalized on quantile (NormalizeBetweenArrays), otherwise no normalization will be performed.\n")
        }
    } else {
        stop("\n- Option technology is mandatory and specifies data type (technology). Options are; array, seq, ms or other. See user manual for specifics.\n")
    }
    return(data)
}

