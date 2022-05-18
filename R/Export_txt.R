#' @title Export DEA results into .txt
#' @description A function for exporting DEA results into a .txt file.
#' @param res.DEA A list of data frames from DEA_feature_apply containing gene counts, p-values, FDRs and logFCs.
#' @param filename The name of the output file.
#' @export
#' @import
#' @seealso
#' @return .txt table
#' @examples \dontrun{
#' ...
#' }

TextOutput <- function(res.DEA, filename) {
    if (is.null(res.DEA)) {
        cat("\nDifferential Expression/Abundance Analysis yielded no results. Is your logFC or FDR cut-offs too strict?\n")
    } else {
        res.DEA <- do.call(rbind, unlist(res.DE, recursive=FALSE))
        my.names <- gsub("1[.](.*)|2[.](.*)", "", rownames(res.DEA))
        my.names <- gsub("group1|group2", "", my.names)
        res.DEA$comparison <- my.names
        file <- try(write.table(res.DEA, paste0(filename,".txt"), sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE), silent = TRUE)
        if (class(file) == "try-error") {
            stop("\n- DEA yielded no results. Aren't your logFC or FDR cut-offs too strict?")
        }
    }
    return(res.DEA)
}
