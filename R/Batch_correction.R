#' @title Batch correction
#' @description Function for batch correction
#' @param data a gene counts matrix
#' @param databatch TRUE/FALSE value indicating activation/deactivation of batch correction. This is automatically set based on the definition of the "batches" parameter.
#' @param batch a factor derived from metadata column including information about a batch for each sample from data
#' @param group a data.frame derived from metadata column selected as a sample group (eg. diagnosis)
#' @param technology a string vector of length 1 (or two in case of 2 datasets) defining technology used for generating the data. Allowed types are: 'array', 'seq', 'ms' or 'other'.
#' @import sva
#' @export
#' @seealso
#' @return a data.frame of batch corrected gene counts
#' @examples \dontrun{
#' ...
#' }

batchCorrect<-function(data,batch,databatch=FALSE,group,technology){

    data.batch=NULL
    if (databatch == TRUE){
        if (length(batch) > 0) {
            design <-  model.matrix(~group)

            if (technology == "seq") {
                data.batch <- ComBat(as.matrix(data$E), batch, design, par.prior=TRUE,prior.plots=FALSE)
            } else {
                data.batch <- ComBat(as.matrix(data), batch, design, par.prior=TRUE,prior.plots=FALSE)
            }

        }
        if (is.null(batch)) {
            print("Batch correction selected but 'batch' parameter remains NULL ")
        }
        if (ncol(data)!=length(batch)) {
            print("Batch correction selected but numbers of samples and provided batch metadata don't fit")
        }
    }

    return(list("data.batch"=data.batch))
}
