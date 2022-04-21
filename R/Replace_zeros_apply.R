#' @title Check of negative values and zeros
#' @description Checking the presence of negative and zero values in the data and replacing zeros.
#' @param data a dataframe of expression/abundance counts
#' @param group a factor specifying sample group from metadata
#' @export
#' @seealso
#' @return a data frame without zero counts
#' @examples \dontrun{
#' ...
#' }

##note which was added in the data parameter description: N.B only a subset of variables should be input, not intended for the full expression matrix! 
##WHY?

ReplaceZerosApply<-function(data=NULL,group=NULL){

    # Check if data contains zeros and negative values.
    hasZero <- unique(as.vector(data == 0))
    hasNeg <- unique(as.vector(data < 0))

    data.original=NULL
    if (TRUE %in% hasNeg) {
        print("Warning: Dataset includes negative values which might cause zeros to be replaced into negative values.")
    }
    if (TRUE %in% hasZero) {
        print("Dataset include zero values which are being replaced.")
        data.original <- data
        data <- ReplaceZero2(data, group)
        print("Zeros were replaced.")
    }


    return(list("data"=data,"data.original"=data.original))

}
