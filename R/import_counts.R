#' @title Read gene counts
#' @description Reading gene count matrix files in txt or xlsx format. Counts must be a numeric value. Check a sample name format after the import.
#' @param my.data a dataframe of expression/abundance counts.
#' @export
#' @import openxlsx
#' @import data.table
#' @import plyr
#' @import scales
#' @examples \dontrun{
#' ...
#' }
#my.data=pocet
importCounts <- function(my.data) {
        file <- try(my.data <- read.table(my.data, header = TRUE, row.names = 1), silent = TRUE)
        if (class(file) == "try-error") {
            cat("\n- Data file is not .txt, trying .xlsx\n")
            file <- try(my.data <- openxlsx::read.xlsx(my.data, colNames = TRUE, rowNames = TRUE), silent = TRUE)
            if (class(file) == "try-error") {
                stop("\n- Data file must be .txt or .xlsx\n")
            }
        }
        my.data<-as.matrix(my.data)
        # Average duplicates and get IDs
#        colnames(my.data)[1] <- "IDs"
#        IDs <- my.data$IDs
#        my.data$IDs <- NULL

#        my.data <- as.data.frame(lapply(my.data, as.numeric))
#        my.data$IDs <- IDs
#        data.frame(data.table(my.data)[, lapply(.SD, mean), by=IDs])

#        my.names <- my.data$IDs
#        my.data$IDs <- NULL
#        my.data <- as.matrix(my.data)
#        rownames(my.data) <- gsub("-", "_", my.names)

        rm(file)
        return(my.data)
}

