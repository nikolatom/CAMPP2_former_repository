#' @title Read metadata
#' @description Reading metadata file in txt or xlsx format
#' @param my.data a dataframe of metadata table.
#' @export
#' @import openxlsx
#' @import data.table
#' @import plyr
#' @import scales
#' @examples \dontrun{
#' ...
#' }

importMetadata <- function(my.data) {
        file <- try(my.data <- read.table(my.data, header = TRUE), silent = TRUE)
        if (class(file) == "try-error") {
            cat("\n- Metadata file is not .txt, trying .xlsx\n")
            file <- try(my.data <- openxlsx::read.xlsx(my.data, colNames = TRUE, rowNames = FALSE), silent = TRUE)
            if (class(file) == "try-error") {
                stop("\n- Metadata file must be .txt or .xlsx\n")
            }
        } else {
        "invalid data type"
        }
        rm(file)
        return(my.data)
}

