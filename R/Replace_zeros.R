#' @title ReplaceZeros
#' @description Checking Data for NA Values.
#' @param my.data a dataframe of expression/abundance counts.
#' @param my.group a vector of integers specifying group; should be represented by a column in metadata file.
#' @export
#' @import impute
#' @seealso
#' @return Dataframe with replaced zeros
#' @examples \dontrun{
#' ReplaceZero (counts, metadata[1])
#' }

# #LOADING LIB ONLY FOR TESTING PURPOSES
# library("impute")
# counts2<-ReplaceZero(counts, metadata[,10])
#
# #NOTES: double check if the % of the missing values are OK
# ##TO BE TESTED IN THE LATER STEPS


ReplaceZero <- function(my.data, my.group) {

    smallestGr <- min(as.numeric(table(my.group)))

    greaterthanBG <- apply(my.data, 1, function(x) sum(x > 0))
    lessthanBG  <- which(as.numeric(greaterthanBG) < smallestGr)

    if (length(lessthanBG) > 0) {
        my.data <- my.data[-lessthanBG,]
    }

    min_per_row <- as.vector(apply(my.data, 1, function(x) min(x[x != 0])))
    for(i in 1:nrow(my.data)){
        my.data[i, my.data[i,] == 0] <- min_per_row[i]
    }
    return(my.data)
}
