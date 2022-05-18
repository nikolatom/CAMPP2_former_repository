#' @title Text Output Function
#' @description Function to export results from differential expression analysis into a .txt file.
#' @param res.DE A list of data frames from DA_feature_apply containing gene counts, p-values, FDRs and logFC.
#' @param filename The name of the output file.
#' @export
#' @import
#' @seealso
#' @return .txt table
#' @examples \dontrun{
#' ...
#' }

TextOutput <- function(res.DE, filename) {
    if (is.null(res.DE)) {
        cat("\nDifferential Expression/Abundance Analysis yielded no results. Is your logFC or FDR cut-offs too strict?\n")
    } else {
        res.DE <- do.call(rbind, unlist(res.DE, recursive=FALSE))
        my.names <- gsub("1[.](.*)|2[.](.*)", "", rownames(res.DE))
        my.names <- gsub("arg.group|arg.sgroup", "", my.names)
        res.DE$comparison <- my.names
        file <- try(write.table(res.DE, paste0(filename,".txt"), sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE), silent = TRUE)
        if (class(file) == "try-error") {
            stop("\n- Differential Expression/Abundance Analysis yielded no results. Is your logFC or FDR cut-offs too strict? Also, check you metadata file has the correct minimum required columns for analysis.")
        }
    }
    return(res.DE)
}
