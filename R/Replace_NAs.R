#' @title Replace "NA" values
#' @description The features (rows) and samples (columns) with high proportion of "NAs" will be removed based on defined thresholds (70% for rows and 80% for columns by default). Remaining "NA" values will be replaced by imputed values. Missing value imputation is done using llsImpute (1st round). If 1st round wasn't successful or resulted in negative values, impute.knn algorithm is applied (rows with more than 50% missing values will be imputed using the overall mean per sample).
#' @param data a dataframe of gene/abundance counts.
#' @param pct.NA.row a number defining maximal percentage of NA values per row (feature). Rows with a higher percentage of NA values will be removed (70 by default).
#' @param pct.NA.column a number defining maximal percentage of NA values per column (sample). Columns with a higher percentage of NA values will be removed (80 by default).
#' @export
#' @import impute
#' @import pcaMethods
#' @seealso
#' @return Data frame with replaced "NA" values
#' @examples \dontrun{
#' ...
#' }

ReplaceNAs <- function(data,pct.NA.row=70,pct.NA.column=80) {
    na_row <- apply(data, 1, function(x) (sum(is.na(x))/ncol(data))*100)
    cat(paste0(" The input data has between " , round(min(na_row), digits = 2), "% - ", round(max(na_row), digits = 2),"%", " missing (NA) values per row. Features (rows) with more than ",pct.NA.row,"% missing values will be removed. \n"))
    removeNA=NULL
    removeNA <- which(as.vector(na_row) > pct.NA.row)
    if(length(removeNA)>0){
        cat(paste0(length(removeNA),"lines will be removed because of a high percentage of NAs"))
        data <- data[-removeNA,]
    } else {
        print("All rows fill requirement on percentage of NAs")
    }

    na_col <- apply(data, 2, function(x) (sum(is.na(x))/nrow(data))*100)
    cat(paste0(" The input data has between " , round(min(na_col), digits = 2), "% - ", round(max(na_col), digits = 2),"%", " missing (NA) values per column. Samples (columns) with more than ", pct.NA.column,"% missing values will be removed. \n N.B. Uncertainty increases with number of missing values! \n "))
    removeNA=NULL
    removeNA <- which(as.vector(na_col) > pct.NA.column)
    if(length(removeNA)>0){
        cat(paste0(length(removeNA),"columns will be removed because of a high percentage of NAs"))
        data <- data[,-removeNA]
    } else {
        print("All columns fill requirement of maximum percentage of NAs in the rows and columns.")
    }

    still.NA <- c(unique(as.vector(is.na(data))))
    if (TRUE %in% still.NA) {
        print("Dataset still contains NA values which will be replaced.")
        varnames <- rownames(data)

        if (checkData(as.matrix(data))[1] == FALSE) {
            data <- as.data.frame(lapply(data, as.numeric))
        }

        data.lls=NULL
        do_knn=FALSE

        print("Running Missing value estimation using local least squares (llsImpute); k=10, correlation=spearman.")

        file <- try(data.lls <- data.frame(completeObs(llsImpute(as.matrix(data), k = 10, correlation="spearman", allVariables=TRUE))), silent =TRUE)

        if (class(file) != "try-error") {
            hasNegB <- unique(as.vector(data < 0))
            hasNegA <- unique(as.vector(data.lls < 0))
            do_knn <- TRUE %in% hasNegA & FALSE %in% hasNegB
        } else {
            do_knn <- TRUE
        }

        if (do_knn==TRUE){
            print("The first round of missing values imputation failed or resulted in negative values. Running additional missing value imputation using impute.knn. Rows with more than 50% missing values will be imputed using the overall mean per sample.")
            data <- impute.knn(as.matrix(data), rowmax = 0.5)                           ##This might generate "Error: C stack usage  xxx is too close to the limit"
            data <- data.frame(data$data)
        } else {
            print("Missing value estimation using local least squares finished successfully.")
            data <- data.lls
        }
        rownames(data) <- varnames
    }

    return(data)
}


