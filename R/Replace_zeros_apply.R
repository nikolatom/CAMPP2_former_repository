#' @title Apply detection and replacement zero values in the data
#' @description Checking the presence of zero values in the data and replacing them.
#' @param data1 a dataframe of expression/abundance counts, N.B only a subset of variables should be input, not intended for the full expression matrix!
#' @param data2 a dataframe of expression/abundance counts for a second dataset
#' @param group1 a factor specifying sample group from metadata1
#' @param group2 a factor specifying sample group from metadata2
#' @export
#' @seealso
#' @return a data frame without zero counts
#' @examples \dontrun{
#' ...
#' }


ReplaceZerosApply<-function(data1,data2=NULL,group1,group2=NULL){

    # Check if data contains zeros and negative values.
    hasZeroD <- unique(as.vector(data1 == 0))
    hasNegD <- unique(as.vector(data1 < 0))

    data1.original=NULL
    if (TRUE %in% hasNegD) {
        print("data1 includes negative values which causes zeros to become negative values")
    }
    if (TRUE %in% hasZeroD) {
        print("data1 include zero values which will be replaced")
        data1.original <- data1
        data1 <- ReplaceZero(data1, group1)
    }


    hasZeroS=NULL
    if (!is.null(data2)){
        hasZeroS <- unique(as.vector(data2 == 0))
        hasNegS <- unique(as.vector(data2 < 0))
    }

    data2.original=NULL
    if (TRUE %in% hasNegS) {
        print("data2 includes negative values which causes zeros to become negative values")
    }
    if (TRUE %in% hasZeroS) {
        print("data2 include zero values which will be replaced")
        data2.original <- data2
        data2 <- ReplaceZero(data2, group2)
    }

    return(list("data1"=data1,"data2"=data2,"data1.original"=data1.original,"data2.original"=data2.original))

}
