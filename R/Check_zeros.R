#' @title Detect and replace zero values in the data
#' @description Checking the presence of zero values in the data and replacing them.
#' @param my.data a dataframe of expression/abundance counts, N.B only a subset of variables should be input, not intended for the full expression matrix!
#' @export
#' @seealso
#' @return normalized counts
#' @examples \dontrun{
#' ...
#' }


checkZeros<-function(data,sdata=NULL){

    # Check if data contains zeros and negative values.
    hasZeroD <- unique(as.vector(data == 0))
    hasNegD <- unique(as.vector(data < 0))

    if(transform[1] %in% c("log2", "log10", "logit")) {
        if (TRUE %in% hasNegD) {
            stop("\n- Data contains negative values and cannot be log transformed. Re-run command WITHOUT argument transform  or alternatively if using two datasets, specify 'none' as the transforminput for the dataset with negative values, e.g. 'none,log2' or 'log2,none'.\n")
        } else {
            if (TRUE %in% hasZeroD) {
                data.original <- data
                data <- ReplaceZero(data, group)
            }
        }
    }




    if (!is.null(sdata)){
        hasZeroS <- unique(as.vector(sdata == 0))
        hasNegS <- unique(as.vector(sdata < 0))
    }

    if(!is.null(sdata) & transform[2] %in% c("log2", "log10", "logit")) {
        if (TRUE %in% hasNegS) {
            stop("\n- Second dataset contains negative values and cannot be log transformed. Re-run command WITHOUT argument transform  or alternatively if using two datasets, specify 'none' as the transforminput for the dataset with negative values, e.g. 'none,log2' or 'log2,none'.\n")
        } else {
            if (TRUE %in% hasZeroS) {
                sdata.original <- sdata
                sdata <- ReplaceZero(sdata, group)
            }
        }
    }
    my_list<-list("data"=data,"sdata"=sdata,"data.original"=data.original,"sdata.original"=sdata.original)
    return(my_list)

}
