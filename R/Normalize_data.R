#' @title Normalization and Transformation
#' @description This function normalizes and transforms feature counts depending
#' on data types and selected normalization/standardization methods.
#' Transformation methods include options log2, logit and log10. If technology
#' is "seq", a Voom transform is applied automatically. Standardization methods
#' options include "mean" or "median". If technology is "seq" or "array", this
#' option is ignored, instead, sequencing data are normalized by "TMM";
#' micro-arrays data are normalized by "quantile" standardization.
#' @param data a data frame of expression/abundance counts
#' @param group a factor specifying group for each sample (e.g. could be
#' represented by a column from a metadata file).
#' @param transform a string vector of length 1 defining transformation type for
#' each dataset ("log2", "logit", "log10"). If technology is "seq", this option
#' is ignored an Voom transform is applied..
#' @param standardize a string vector of length 1 defining standardization
#' method ("mean" or "median"). If technology is "seq" or "array", this option
#' is ignored, instead, seq data are normalized by "TMM"; array data are
#' normalized by "quantile" standardization.
#' @param technology a string vector of length 1 defining technology used for
#' generating the data. Allowed types are: 'array', 'seq', 'ms' or 'other'.
#' @export
#' @import fitdistrplus
#' @return Elist (seq technology) or array (array, ms, other technologies) of
#' normalized and transformed feature counts data
#' @examples {
#' ###In this example, data with fixed zeros are used as an input.
#' campp2_brca_1_normalized<-NormalizeData(data=campp2_brca_1_zeroFix,
#' group=campp2_brca_1_meta$diagnosis, standardize="TMM", transform="voom",
#' technology="seq")
#' }
#'

NormalizeData <- function(data,group,transform,standardize,technology){

    # section in which we define correctness of input parameters
    if (!(technology) %in% c("seq", "array", "ms", "other")) {
        stop("Defined technology is not supported.")
    }

    if (technology != "seq") {
        if (!(transform) %in% c("log2", "logit", "log10")) {
            stop("Defined transformation is not supported")
    	}
    }

    if (!(technology) %in% c("seq", "array")){
        if ( !(standardize) %in% c("mean", "median")) {
            stop("Defined normalization is not supported.")
        }
    }

    #check for negative values
    hasNeg <- unique(as.vector(data < 0))
    hasZero <- unique(as.vector(data == 0))
    if(transform %in% c("log2", "log10", "logit")) {
        if (TRUE %in% hasNeg || TRUE %in% hasZero) {
            stop("\n- Data contains zero/negative values and cannot be log transformed. \n")
        }
    }

    # for RNAseq, perform Voom transform and normalize by TMM
    if (technology == "seq") {
        cat("\n Data are being filtered for lowly expressed variables, normalized using TMM and voom transformed without any other option for normalization and transformation. Options transform and standardize are being ignored.\n")
        data <- DGEList(counts=data)
        design <- model.matrix(~0+group)
        keep <- filterByExpr(data, design)
        data <- data[keep,,keep.lib.sizes=FALSE]
        data <- calcNormFactors(data, method = "TMM")
        data <- voom(data, design, plot=FALSE)

	return(data)
    # this else is not strictly necessary but good for clarity
    } else {

    ##Should we support also another normalization types for "seq" technology?
        if (transform == "log2") {
            cat("\n Data are being log2 transformed.\n")
            data <- log2(data)
        } else if (transform == "logit") {
            cat("\n Data are being logit transformed.\n")
            data <- logit(data)
        } else if (transform == "log10"){
            cat("\n Data are being log10 transformed.\n")
            data <- log10(data)
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
