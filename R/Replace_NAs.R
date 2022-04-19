#' @title Replace "NA" values
#' @description Replacing "NA" values in the data (counts).
#' @param data a dataframe of gene/abundance counts.
#' @export
#' @import impute
#' @import pcaMethods
#' @seealso
#' @return Data frame with replaced "NA" values
#' @examples \dontrun{
#' ...
#' }

ReplaceNAs <- function(data) {
    na_row <- apply(data, 1, function(x) (sum(is.na(x))/ncol(data))*100)
    cat(paste0(" The input data has between " , round(min(na_row), digits = 2), "% - ", round(max(na_row), digits = 2),"%", " missing values per row. Samples (rows) with more than 70% missing values will be removed. \n"))
    removeNA=NULL
    removeNA <- which(as.vector(na_row) > 70)
    if(length(removeNA)>0){
        data <- data[-removeNA,]
    }
    
    na_col <- apply(data, 2, function(x) (sum(is.na(x))/nrow(data))*100)
    cat(paste0(" The input data has between " , round(min(na_col), digits = 2), "% - ", round(max(na_col), digits = 2),"%", " missing values per column. Samples (rows) with more than 80% missing values will be removed. \n N.B. Uncertainty increases with number of missing values! \n "))
    removeNA=NULL
    removeNA <- which(as.vector(na_col) > 80)
    if(length(removeNA)>0){
        data <- data[,-removeNA]
    }
    
    still.NA <- c(unique(as.vector(is.na(data))))
    if (TRUE %in% still.NA) {
        varnames <- rownames(data)
        ###NOT SURE ABOUT THIS PART, I think it wasn't working correctly. I have added 1 line below
        # if (checkData(as.matrix(my.data))[1] == FALSE) {
        #     my.data <- as.data.frame(lapply(my.data, as.numeric))
        # }
        
        data<-as.data.frame(lapply(data, as.numeric))
        data.lls=NULL
        print("Running Missing value estimation using local least squares (llsImpute)")
        file <- try(data.lls <- data.frame(completeObs(llsImpute(as.matrix(data), k = 10, correlation="spearman", allVariables=TRUE))), silent =TRUE)
        
        hasNegB <- unique(as.vector(data < 0))
        hasNegA <- unique(as.vector(data.lls < 0))
        if (class(file) == "try-error" || TRUE %in% hasNegA & hasNegB == FALSE) {       ##If the 1st approach not working, try impute.knn)
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


