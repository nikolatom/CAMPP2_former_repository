#' @title Replace zeros
#' @description Replacing zero values in the data (counts).
#' @param data a dataframe of gene/abundance counts.
#' @param group a factor specifying group which should be represented by a column from a metadata file.
#' @export
#' @import impute
#' @seealso
#' @return a data frame with replaced zeros
#' @examples \dontrun{
#' ...
#' }


ReplaceZero <- function(data, group) {

    smallestGr <- min(as.numeric(table(group))) # a size of the smallest group of samples
    greaterthanBG <- apply(data, 1, function(x) sum(x > 0)) #rows with sum>0
    lessthanBG  <- which(as.numeric(greaterthanBG) < smallestGr)

    if (length(lessthanBG) > 0) {
        data <- data[-lessthanBG,] #removal of lines with low counts
    }

    min_per_row <- as.vector(apply(data, 1, function(x) min(x[x != 0])))
    for(i in 1:nrow(data)){
        data[i, data[i,] == 0] <- min_per_row[i] #substitution of zeros with min values per row
    }
    return(data)
}
