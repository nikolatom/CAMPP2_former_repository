#' @title Replace zeros
#' @description A function for a detection of zero/negative values and for replacing zero values in the data.
#' @param data a dataframe of gene/abundance counts.
#' @param group a factor specifying group which should be represented by a column from a metadata file.
#' @export
#' @import impute
#' @seealso
#' @return a data frame with replaced zeros and a dataframe with original data
#' @examples \dontrun{
#' ...
#' }


ReplaceZeros <- function(data, group) {
    hasZero <- unique(as.vector(data == 0))
    if(TRUE %in% hasZero){
        print("data includes 0-value(s)")
    }else{
        print("data doesn't include 0-value(s)")
    }

    hasNeg <- unique(as.vector(data < 0))
    if(TRUE %in% hasNeg){
        print("data includes negative value(s)")
    }else{
        print("data doesn't include negative value(s)")
    }


    smallestGr <- min(as.numeric(table(group))) # a size of the smallest group of samples
    greaterthanBG <- apply(data, 1, function(x) sum(x > 0)) #counting non-zero counts for each feature
    lessthanBG  <- which(as.numeric(greaterthanBG) < smallestGr) #features having number of zero counts higher than the size of the smallest sample group.

    data.original<-NULL
    if (length(lessthanBG) > 0) {
        data.original<-data
        data <- data[-lessthanBG,] #removal of lines with low counts
        print(paste0(length(lessthanBG)," features will be removed because of the low gene counts"))
    }

    min_per_row <- as.vector(apply(data, 1, function(x) min(x[x != 0])))
    for(i in 1:nrow(data)){
        data[i, data[i,] == 0] <- min_per_row[i] #substitution of zeros with min values per row
    }
    return(list(data, data.original))
}
