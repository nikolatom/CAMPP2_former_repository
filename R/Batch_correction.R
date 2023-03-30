#' @title Batch correction
#' @description This function removes batch effects (using ComBat function) which could be causing a significant heterogeneity across batches of data. Batch corrected data are intended for explanatory purposes, NOT for a DEA analysis.
#' @param data an Elist ("seq" technology) or a matrix ("array", "ms", "other" technologies) of normalized and transformed feature (gene) counts data
#' @param batch a factor specifying batch for each sample (e.g. could be represented by a column from a metadata file)
#' @param group a factor specifying group for each sample (e.g. could be represented by a column from a metadata file)
#' @param technology a string vector of length 1 defining technology used for generating the data. Allowed types are: "array", "seq", "ms" or "other"
#' @import sva
#' @export
#' @return a list including matrix (array) of batch corrected feature counts
#' @examples {
#' campp2_brca_1_batchCorrected<-BatchCorrect(data=campp2_brca_1_normalized,
#' batch=campp2_brca_1_meta$tumor_stage,group=campp2_brca_1_meta$diagnosis,
#' technology="seq")
#' }

BatchCorrect<-function(data,batch,group,technology){

    if (!(technology) %in% c("seq", "array", "ms", "other")) {
        stop("Defined technology is not supported.")
    }

    data.batch=NULL
    if (ncol(data)!=length(batch)) {
        stop("Batch correction was selected BUT a number of samples and length of provided batch metadata don't fit")
    }
    else if (length(batch) > 0) {
        design <-  model.matrix(~group)
        if (technology == "seq") {
            data.batch <- ComBat(as.matrix(data$E), batch, design, par.prior=TRUE,prior.plots=FALSE)
        } else {
            data.batch <- ComBat(as.matrix(data), batch, design, par.prior=TRUE,prior.plots=FALSE)
            }
    }

    return(data.batch)
}
