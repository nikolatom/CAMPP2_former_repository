#' @title Fix zeros
#' @description This function detects zero/negative values and replaces zero
#' values in the data. In the first step, the function checks for the presence
#' of zero and negative values in the data. In the second step, (activated by
#' default), features having sum of zero counts higher than the size of the
#' smallest sample group will be removed (e.g. genes represented by zero
#' counts/transcripts in more than n samples will be removed; n = size of the
#' smallest group). Next, remaining zeros are substituted with minimal values
#' (>0) observed in each feature across all samples.
#' @param data a dataframe of gene/abundance counts.
#' @param group group a factor specifying samples' group (e.g. could be
#' represented by a column from a metadata file).
#' @param remove.sparse.features a logical argument (TRUE/FALSE) for removal of
#' features with sum of zero counts larger than the size of the smallest
#' sample group, and for replacing the remaining zeros.  Default is TRUE.
#' @export
#' @import impute
#' @return a data frame with fixed zeros. Features having sum of zero counts
#' higher than the size of the smallest sample group are removed and
#' remaining zeros will be replaced by default.
#' @examples {
#' ###In this example, data with fixed NA values are used as an input.
#' FixZeros(data=campp2_brca_1_replacedNAs,group=campp2_brca_1_meta$diagnosis, remove.sparse.features=TRUE)
#' }


FixZeros <- function(data, group, remove.sparse.features=TRUE) {
    if(any(data==0)){
        print("data includes 0-value(s)")
    }else{
        print("data doesn't include 0-value(s)")
    }

    if(any(data < 0)){
        print("data includes negative value(s)")
    }else{
        print("data doesn't include negative value(s)")
    }

    ###Removal of features with high 0-counts
    if(remove.sparse.features==TRUE){
        smallestGr <- min(as.numeric(table(group))) # a size of the smallest
        #group of samples
        greaterthanBG <- apply(data, 1, function(x) sum(x > 0)) #counting
        #non-zero counts for each feature
        lessthanBG  <- which(as.numeric(greaterthanBG) < smallestGr) #features
        #having number of zero counts higher than the size of the smallest
        #sample group.

        if (length(lessthanBG) > 0) {
            data <- data[-lessthanBG,] #removal of lines with low counts
            print(paste0(length(lessthanBG)," features will be removed because of the low counts"))
        }
    }

    ###Replacement of zeros by the lowest value for each feature
    min_per_row <- as.vector(apply(data, 1, function(x) min(x[x != 0])))
    for(i in 1:nrow(data)){
        data[i, data[i,] == 0] <- min_per_row[i] #substitution of zeros with
        #min values per row
    }

    return(data)
}
