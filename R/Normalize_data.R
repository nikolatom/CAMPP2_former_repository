#' @title Normalization
#' @description Normalizing features' counts depending on data types and selected normalization/standardization methods
#' @param my.data a dataframe of expression/abundance counts, N.B only a subset of variables should be input, not intended for the full expression matrix!
#' @param my.technology Possible values are: "array" (mircoarray data), "seq" (high throughput sequencing data), "ms" (mass spectrometry data) or "other" (other type).
#' @param my.group a vector of integers specifying group; should be represented by a column in metadata file.
#' @param my.transform a normalization type - "log2", "logit" or "log10"
#' @param my.standardize a standardization method - "mean" or "median"
#' @param my.data.original = NULL
#' @export
#' @import fitdistrplus
#' @seealso LINK TO MANUAL (TRANSFORMATION/NORMALIZATION)
#' @return normalized counts
#' @examples \dontrun{
#' ...
#' }

# my.technology=technology
# my.transform=transform
# my.data=dataset1
# my.standardize=standardize
# my.group=group
NormalizeData <- function(my.technology, my.data, my.group, my.transform, my.standardize, my.data.original = NULL) {
    if (my.technology == "seq") {
        if (!is.null(my.data.original)) {
            my.data <- my.data.original
        }
        my.data <- DGEList(counts=my.data)
        design <- model.matrix(~0+my.group)
        keep <- filterByExpr(my.data, design)
        my.data <- my.data[keep,,keep.lib.sizes=FALSE]
        my.data <- calcNormFactors(my.data, method = "TMM")
        my.data <- voom(my.data, design, plot=TRUE)
        cat("\n Data will be filtered for lowly expressed variables, normalized and voom transformed.\n")

    } else if (my.technology %in% c("array", "ms", "other")) {
        if (my.transform == "log2") {
            my.data <- log2(my.data)
        } else if (my.transform == "logit") {
            my.data <- logit(my.data)
        } else if (my.transform == "log10"){
            my.data <- log10(my.data)
        } else {
            cat("\n my.transform is not specified for data, log transformation will NOT be performed.\n")
            my.data <- my.data
        }
        if (my.standardize == "mean") {
            my.data <- scale(my.data, scale = FALSE)
            NB <- " N.B This pipeline does not handle background correction of single-channel intensity data or within array normalization two-color intensity data. See limma manual section on array normalization for more on this. Data may be fully normalized with limma (R/Rstudio) or another software and the pipeline re-run."
            cat(paste0("\n my.technology = array and my.standardize = mean. Data will be mean centered.", NB, "\n"))

        } else if (my.standardize == "median") {
            rowmed <- apply(my.data,1,median)
            my.data <- my.data - rowmed
            NB <- " N.B This pipeline does not handle background correction of single-channel intensity data or within array normalization two-color intensity data. See limma manual section on array normalization for more on this. Data may be fully normalized with limma (R/Rstudio) or another software and the pipeline re-run."
            cat(paste0("\n my.technology = array and my.standardize = median. Data will be median centered.", NB, "\n"))
        } else if (!(my.standardize %in% c("mean", "median")) & my.technology == "array") {
            my.data <- normalizeBetweenArrays(my.data, method="quantile")
            NB <- " N.B This pipeline does not handle background correction of single-channel intensity data or within array normalization two-color intensity data. See limma manual section on array normalization for more on this. Data may be fully normalized with limma (R/Rstudio) or another software and the pipeline re-run."
            cat(paste0("\n my.technology = array. Data will be normalized on the quantiles.",NB, "\n"))
        } else {
            cat("\n- No standardization requested. If argument my.technology is 'array', data will be normalized on quantile (NormalizeBetweenArrays), otherwise no normalization will be performed.\n")
        }
    } else {
        stop("\n- Option my.technology is mandatory and specifies data type (technology). Options are; array, seq, ms or other. See user manual for specifics.\n")
    }
    return(my.data)
}
