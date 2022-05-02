#' @title Normalization and Transformation
#' @description This function normalizes and transforms feature counts depending on data types and selected normalization/standardization methods.
#' @param data a dataframe of expression/abundance counts
#' @param group a factor specifying sample group from metadata
#' @param data.original a original data (before removal of zeros)
#' @param transform a string vector of length 1 defining transformation type for each dataset ("log2", "logit", "log10", NULL for no transformation). If technology is "seq", this option is ignored an Voom transform is applied instead.
#' @param standardize a string vector of length 1 defining standardization method  ("mean" or "median"). If technology is "seq" or "array", this option is ignored, instead, seq data is normalized by "TMM"; array data is normalized by "quantile" standardization
#' @param technology a string vector of length 1 (or two in case of 2 datasets) defining technology used for generating the data. Allowed types are: 'array', 'seq', 'ms' or 'other'.
#' @export
#' @import fitdistrplus
#' @seealso
#' @return Elist (seq technology) or array (array, ms, other technologies) of normalized, filtered and transformed (gene) counts data
#' @examples \dontrun{
#' ...
#' }
#'

NormalizeData <- function(data,data.original,group,transform,standardize,technology) {
 
    # section in which we define correctness of input parameters
    if (!(technology) %in% c("seq", "array", "ms", "other")) {
        stop("Defined technology is not supported.")
    }

    if ((technology != "seq") {
        if (!(transform) %in% c("log2", "logit", "log10", NULL)) {
            stop("Defined transformation is not supported")
	}
    }

    if (!(technology) %in% c("seq", "array")){
        if ( !(standardize) %in% c("mean", "median", NULL)) {
            stop("Defined normalization is not supported.")
        }
    }



    # for RNAseq, perform Voom transform and normalize by TMM
    if (technology == "seq") {
        if (!is.null(data.original)) {
            data <- data.original #takes the original data (before removing zeros)
        }
        cat("\n Data are being filtered for lowly expressed variables, normalized using TMM and voom transformed without any other option for normalization and transformation. Options transform and standardize are being ignored.\n")
        data <- DGEList(counts=data)
        design <- model.matrix(~0+group)
        keep <- filterByExpr(data, design)
        data <- data[keep,,keep.lib.sizes=FALSE]
        data <- calcNormFactors(data, method = "TMM")
        data <- voom(data, design, plot=TRUE)

	return(data)	
    # this else is not strictly necessary but good for clarity
    } else {

    ##Should we support also another normalization types for "seq" technology
        if (transform == "log2") {
            cat("\n Data are being log2 transformed.\n")
            data <- log2(data)
        } else if (transform == "logit") {
            cat("\n Data are being logit transformed.\n")
            data <- logit(data)
        } else if (transform == "log10"){
            cat("\n Data are being log10 transformed.\n")
            data <- log10(data)
        } else if (is.null(transform)){
            cat("\n Transformation will NOT be performed.\n")
        }
    }
        
    if (technology == "array") {
        cat("\n Data will be normalized on the quantiles (default for microarray technology). Option standardize will be ignored.\n")
        data <- normalizeBetweenArrays(data, method="quantile")
    } else {
        if (standardize == "mean") {
            cat("\n Data will be mean centered. \n")
            data <- scale(data, scale = FALSE)
        } else if (standardize == "median") {
            cat("\n Data will be median centered. \n")
            rowmed <- apply(data,1,median)
            data <- data - rowmed
        }
    }
    return(data)
}
