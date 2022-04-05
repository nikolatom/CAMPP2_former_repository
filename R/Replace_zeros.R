#' @title Replace zeros
#' @description Replaces zero values in the data (counts).
#' @param data a dataframe of expression/abundance counts.
#' @param group a factor specifying group; should be represented by a column in metadata file.
#' @export
#' @import impute
#' @seealso
#' @return Data frame with replaced zeros
#' @examples \dontrun{
#' ...
#' }


ReplaceZero <- function(data, group) {

    smallestGr <- min(as.numeric(table(group)))

    greaterthanBG <- apply(data, 1, function(x) sum(x > 0))
    lessthanBG  <- which(as.numeric(greaterthanBG) < smallestGr)

    if (length(lessthanBG) > 0) {
        data <- data[-lessthanBG,]
    }

    min_per_row <- as.vector(apply(data, 1, function(x) min(x[x != 0])))
    for(i in 1:nrow(data)){
        data[i, data[i,] == 0] <- min_per_row[i]
    }
    return(data)
}
