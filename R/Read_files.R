#' @title Read files
#' @description Reading input files
#' @param my.expr This param is set to TRUE or FALSE by CAMPP2
#' @param my.data a dataframe of expression/abundance counts.
#' @export
#' @import openxlsx
#' @import data.table
#' @import plyr
#' @import scales
#' @examples \dontrun{
#' ...
#' }

ReadMyFile <- function(my.data, my.expr) {
    if(my.expr == TRUE) {
        file <- try(my.data <- openxlsx::read.xlsx(my.data, colNames = TRUE, rowNames = FALSE), silent = TRUE)
        if (class(file) == "try-error") {
            cat("\n- Data file is not .xlsx, trying .txt\n")
            file <- try(my.data <- read.delim(my.data, header = TRUE), silent = TRUE)
            if (class(file) == "try-error") {
                file <- try(my.data <- read.delim(my.data, header = TRUE, sep=";"), silent = TRUE)
                if (class(file) == "try-error") {
                    stop("\n- Data file must be .xlsx or .txt\n")
                }
            }
        }




        # Average duplicates and get IDs
        colnames(my.data)[1] <- "IDs"
        IDs <- my.data$IDs
        my.data$IDs <- NULL

        my.data <- as.data.frame(lapply(my.data, as.numeric))
        my.data$IDs <- IDs
        my.data <-  data.frame(data.table(my.data)[, lapply(.SD, mean), by=IDs])

        my.names <- my.data$IDs
        my.data$IDs <- NULL
        my.data <- as.matrix(my.data)
        rownames(my.data) <- gsub("-", "_", my.names)

    } else if(data.type == "metadata"){
        file <- try(my.data <- read.xlsx(my.data, colNames = TRUE, rowNames = FALSE), silent = TRUE)
        if (class(file) == "try-error") {
            cat("\n- Metadata file is not .xlsx, trying .txt\n")
            file <- try(my.data <- read.delim(my.data, header = TRUE), silent = TRUE)
            if (class(file) == "try-error") {
                file <- try(my.data <- read.delim(my.data, header = TRUE, sep =";"), silent = TRUE)
                if (class(file) == "try-error"){
                    stop("\n- Metadata file must be .xlsx or .txt\n")
                }
            }
        }
    } else {
        "invalid data type"
    }
    rm(file)
    return(my.data)
}

