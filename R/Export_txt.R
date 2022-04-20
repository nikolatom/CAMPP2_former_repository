#' @title Text Output Function
#' @description Function to export results into a .txt file
#' @param my.list A list of dataframes from DA_feature_apply with p-values, FDRs, logFC ect.
#' @param my.filename The name of output txt file
#' @export
#' @import
#' @seealso
#' @return .txt table
#' @examples \dontrun{
#' ...
#' }

TextOutput <- function(my.list, my.filename) {
    if (is.null(my.list)) {
        cat("\nDifferential Expression/Abundance Analysis yielded no results. Is your logFC or FDR cut-offs too strict?\n")
    } else {
        my.list <- do.call(rbind, unlist(my.list, recursive=FALSE))
        my.names <- gsub("1[.](.*)|2[.](.*)", "", rownames(my.list))
        my.names <- gsub("arg.group|arg.sgroup", "", my.names)
        my.list$comparison <- my.names
        file <- try(write.table(my.list, paste0(my.filename,".txt"), sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE), silent = TRUE)
        if (class(file) == "try-error") {
            stop("\n- Differential Expression/Abundance Analysis yielded no results. Is your logFC or FDR cut-offs too strict? Also, check you metadata file has the correct minimum required columns for analysis.")
        }
    }
    return(my.list)
}

