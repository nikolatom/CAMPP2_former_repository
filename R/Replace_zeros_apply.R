#' @title Detection of negative values and zeros and replacement of zero values in the data
#' @description Checking the presence of negative and zero values in the data and replacing zeros.
#' @param data a dataframe of expression/abundance counts, N.B only a subset of variables should be input, not intended for the full expression matrix!
#' @param group a factor specifying sample group from metadata
#' @export
#' @seealso
#' @return a data frame without zero counts
#' @examples \dontrun{
#' ...
#' }


ReplaceZerosApply<-function(data=NULL,group=NULL){

    # Check if data contains zeros and negative values.
    hasZero <- unique(as.vector(data == 0))
    hasNeg <- unique(as.vector(data < 0))

    data.original=NULL
    if (TRUE %in% hasNeg) {
        print("Dataset includes negative values which might cause zeros to be replaced into negative values")
    }
    if (TRUE %in% hasZero) {
        print("Dataset include zero values which are being replaced")
        data.original <- data
        data <- ReplaceZero(data, group)
    }


    return(list("data"=data,"data.original"=data.original))

}
