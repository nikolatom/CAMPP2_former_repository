#' @title Batch correction
#' @description Function for a batch correction. Batch corrected data are NOT intended for a DEA analysis with limma.
#' @param data Elist (seq technology) or array (array, ms, other technologies) of normalized, filtered and transformed (gene) counts data
#' @param databatch TRUE/FALSE value indicating activation/deactivation of batch correction. Boolean value is automatically set based on the definition of the "batches" parameter.
#' @param batch a factor derived from metadata column including information about a batch for each sample from data
#' @param group a factor derived from metadata column selected as a sample group (e.g. diagnosis)
#' @param technology a string vector of length 1 (or two in case of 2 datasets) defining technology used for generating the data. Allowed types are: 'array', 'seq', 'ms' or 'other'.
#' @import sva
#' @export
#' @seealso
#' @return a list including matrix (array) of batch corrected feature counts
#' @examples \dontrun{
#' ...
#' }

BatchCorrect<-function(data,batch,databatch=FALSE,group,technology){
    
    data.batch=NULL
    if (ncol(data)!=length(batch)) {
        stop("Batch correction selected BUT number of samples and length of provided batch metadata don't fit")
    }
    else if (databatch == TRUE){
        if (length(batch) > 0) {
            design <-  model.matrix(~group)

            if (technology == "seq") {
                data.batch <- ComBat(as.matrix(data$E), batch, design, par.prior=TRUE,prior.plots=FALSE)
            } else {
                data.batch <- ComBat(as.matrix(data), batch, design, par.prior=TRUE,prior.plots=FALSE)
            }
        }
    }

    return(list("data.batch"=data.batch))
}
