#' @title Distribution check
#' @description Function for checking gene counts distribution
#' @param data a gene counts matrix
#' @param sdata a gene counts matrix for a secend dataset
#' @param datacheck TRUE/FALSE parameter for distributional check
#' @param databatch TRUE/FALSE parameter for batch correction
#' @param sdatabatch TRUE/FALSE parameter for batch correction for the second dataset
#' @param data.batch batch corrected gene counts
#' @param sdata.batch batch corrected gene counts for the second dataset
#' @import
#' @export
#' @seealso
#' @return batch corrected gene counts
#' @examples \dontrun{
#' ...
#' }

checkDistr<-function(data,sdata=NULL,datacheck,databatch,sdatabatch,data.batch=NULL,sdata.batch=NULL){
    if (datacheck == TRUE) {
        if (databatch == TRUE) {
            subset.data <- data.batch[sample(nrow(data.batch), 10),]
        } else {
            if (technology[1] == "seq") {
                subset.data <- data$E[sample(nrow(data$E), 10),]
            } else {
                subset.data <- data[sample(nrow(data), 10),]
            }
        }

        list.of.dists <- FitDistributions(subset.data)

        dir.create("DataChecks")
        setwd("DataChecks/")
        PlotDistributions(subset.data, list.of.dists)
        setwd("..")

        rm(subset.data, list.of.dists)
    }



    if (datacheck == TRUE & !is.null(sdata)) {
        if (sdatabatch == TRUE) {
            subset.data <- sdata.batch[sample(nrow(sdata.batch), 10),]
        } else {

            if (technology[2] == "seq") {
                subset.data <- sdata$E[sample(nrow(sdata$E), 10),]
            } else {
                subset.data <- sdata[sample(nrow(sdata), 10),]
            }
        }

        list.of.dists <- FitDistributions(subset.data)

        dir.create("SecondDataChecks")
        setwd("SecondDataChecks/")
        PlotDistributions(subset.data, list.of.dists)
        setwd("..")

        rm(subset.data, list.of.dists)
    }
}
