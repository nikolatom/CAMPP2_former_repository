#' @title ReplaceNAs
#' @description Checking Data for NA Values.
#' @param my.data a dataframe of expression/abundance counts.
#' @export
#' @import impute
#' @seealso
#' @return Dataframe with replaced NA-values
#' @examples \dontrun{
#' ...
#' }

ReplaceNAs <- function(my.data) {
    na_row <- apply(my.data, 1, function(x) (sum(is.na(x))/ncol(my.data))*100)

    cat(paste0("\n- The input data has between " , round(min(na_row), digits = 2), "% - ", round(max(na_row), digits = 2),"%", " missing values per row.\n- Variables (rows) with more than 70% missing values will be removed.\n"))


    removeNA <- which(as.vector(na_row) > 70)
    if(length(removeNA) > 0) {
        my.data <- my.data[-removeNA,]

    }
    na_col <- apply(my.data, 2, function(x) (sum(is.na(x))/nrow(my.data))*100)
    cat(paste0("\n- The input data has between " , round(min(na_col), digits = 2), "% - ", round(max(na_col), digits = 2),"%", " missing values per column.\n- Samples (rows) with more than 80% missing values will be removed. Variables (rows) with more than 50% missing values are imputed using the overall mean per sample.\n... Performing missing value imputation. N.B Uncertainty increases with number of missing values!.\n"))

    removeNA <- which(as.vector(na_col) > 80)
    if(length(removeNA) > 0) {
        my.data <- my.data[,-removeNA]
    }

    still.NA <- unique(as.vector(is.na(my.data)))

    if (TRUE %in% still.NA) {
        varnames <-  rownames(my.data)

        if (checkData(as.matrix(my.data))[1] == FALSE) {
            my.data <- as.data.frame(lapply(my.data, as.numeric))
        }

        file <- try(my.data.lls <- data.frame(completeObs(llsImpute(as.matrix(my.data), k = 10, correlation="spearman", allVariables=TRUE))), silent =TRUE)
        hasNegB <- unique(as.vector(my.data < 0))
        hasNegA <- unique(as.vector(my.data.lls < 0))

        if (class(file) == "try-error" || TRUE %in% hasNegA & hasNegB == FALSE) {
            my.data <- impute.knn(as.matrix(my.data), rowmax = 0.7)
            my.data <- data.frame(my.data$data)
        } else {
            my.data <- my.data.lls
            rm(file)
        }
        rownames(my.data) <- varnames  ###THIS LINE WAS ORIGINALLY BEFORE "return(my.data) but throwing an error - is it correct now?"
    }
    return(my.data)
}


