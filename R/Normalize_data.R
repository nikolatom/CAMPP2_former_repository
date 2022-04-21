#' @title Normalization and Transformation
#' @description This function normalizes and transforms feature counts depending on data types and selected normalization/standardization methods.
#' @param data a dataframe of expression/abundance counts
#' @param group a factor specifying sample group from metadata
#' @param data.original a original data (before removal of zeros)
#' @param transform a string vector of length 1 (or two in case of 2 datasets) defining transformation type for each dataset ("log2", "logit", "log10"). voom transform is automatically applied on "seq" data.
#' @param standardize a string vector of length 1 (or two in case of 2 datasets) defining standardization method  ("mean" or "median"); By default, data from "seq" technology are normalized by "TMM"; "array" technology is standardized using "quantile"
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
  
    if (is.null(technology)){
        stop("\n- Mandatory parameter technology is not defined for a given dataset.\n")
    }
    
    if (technology == "seq") {
        if (!is.null(data.original)) {
            data <- data.original #takes the original data (before removing zeros)
        }
        cat("\n Data are being filtered for lowly expressed variables, normalized using TMM and voom transformed without any other option for normalization and transformation.\n")
        data <- DGEList(counts=data)
        design <- model.matrix(~0+group)
        keep <- filterByExpr(data, design)
        data <- data[keep,,keep.lib.sizes=FALSE]
        data <- calcNormFactors(data, method = "TMM")
        data <- voom(data, design, plot=TRUE)
        }

    ##Should we support also another normalization types for "seq" technology?
    
     else if (technology %in% c("array", "ms", "other")) {
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
            cat("\n Transformation is not specified and will NOT be performed.\n")
            data <- data
        } 
        
        if (standardize == "mean") {
            cat("\n Data will be mean centered test. \n")
            data <- scale(data, scale = FALSE)
        } else if (standardize == "median") {
            cat("\n Data will be median centered. \n")
            rowmed <- apply(data,1,median)
            data <- data - rowmed
        } else if (!(standardize %in% c("mean", "median")) & technology == "array") {
            cat("\n Data will be normalized on the quantiles (default for microarray technology). \n")
            data <- normalizeBetweenArrays(data, method="quantile")
        } else {
            cat("\n Given standardization is not supported or there are no standardization defaults for a given technology.\n")
        } 
         
     } else if (!(technology) %in% c("seq","array","ms","other")){
         stop("Defined technology is not supported.")
     }

    return(data)

}
