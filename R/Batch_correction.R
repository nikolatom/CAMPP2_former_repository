#' @title Batch correction
#' @description Function for batch correction
#' @param data a gene counts matrix
#' @param sdata a gene counts matrix for a secend dataset
#' @param databatch TRUE/FALSE parameter for batch correction
#' @param sdatabatch TRUE/FALSE parameter for batch correction for the second dataset
#' @param batch batch information for each sample
#' @param sbatch batch information for each sample from the second dataset
#' @param group group information for each sample
#' @param sgroup group information for each sample from the second dataset
#' @param technology technology
#' @import sva
#' @export
#' @seealso
#' @return batch corrected gene counts
#' @examples \dontrun{
#' ...
#' }

batchCorrect<-function(data,sdata=NULL,batch=NULL,sbatch=NULL,databatch=NULL,sdatabatch=NULL,group,sgroup=NULL,technology){

    data.batch=NULL
    sdata.batch=NULL
    print("second databatch")
    print(databatch)

    if (databatch == TRUE){
        if (length(batch) > 0) {
            design <-  model.matrix(~group)

            if (technology[1] == "seq") {
                data.batch <- ComBat(as.matrix(data$E), batch, design, par.prior=TRUE,prior.plots=FALSE)
            } else {
                data.batch <- ComBat(as.matrix(data), batch, design, par.prior=TRUE,prior.plots=FALSE)
            }

        } else {
            data.batch <- data
            cat("\n- No column names match specified batches for dataset.\n")
        }
    } else {
        cat("\n- No batch correction requested.\n")
    }


    if (sdatabatch == TRUE){
        if (length(sbatch) > 0) {
            sdesign <- model.matrix(~sgroup)

            if (technology[2] == "seq") {
                sdata.batch <- ComBat(as.matrix(sdata$E), sbatch, sdesign, par.prior=TRUE,prior.plots=FALSE)
            } else {
                sdata.batch <- ComBat(as.matrix(sdata), sbatch, sdesign, par.prior=TRUE,prior.plots=FALSE)
            }
        } else {
            sdata.batch <- sdata
            cat("\n- No column names match specified batches for second dataset. Continuing without batch correction.\n")
        }
    } else {
        cat("\n- No batch correction requested for second dataset.\n")
    }
    my_list<-list("data.batch"=data.batch,"sdata.batch"=sdata.batch)
    return(my_list)
}
