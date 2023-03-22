#' @title Subset DA/DE
#' @description Function which subsets a dataframe of differential expression analysis into a list of dataframes by comparison.
#' @param my.data a dataframe with results of differential expression analysis.
#' @export
#' @return a list of dataframes by comparison
#' @examples {
#' ...
#' }


DFtoList <- function(my.data) {

    my.data$name <- gsub("_", "-", my.data$name)
    comparisons <- levels(as.factor(my.data$comparison))
    df.list <- list()

    for (idx in 1:length(comparisons)) {
        df <- my.data[my.data$comparison == as.character(comparisons[idx]),]
        df.list[[idx]] <- df[,-c(2:4,6)]
    }

    names(df.list) <- comparisons

    df.lengths <- as.numeric(unlist(lapply(df.list, function(x) nrow(x))))
    df.lengths.min <- which(df.lengths < 2)
    if (length(df.lengths.min) > 0) {
        df.list <- df.list[-df.lengths.min]
    }
    return(df.list)
}

