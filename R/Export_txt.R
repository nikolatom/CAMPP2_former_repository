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

ExportDEA <- function(res.DEA, filename) {
    if (is.null(res.DEA)) {
        cat("\nDEA yielded no results. Aren't your logFC or FDR cut-offs too strict?\n")
    } else {
        res.DEA <- do.call(rbind, unlist(res.DEA, recursive=FALSE))
        group.names <- gsub("1[.](.*)|2[.](.*)", "", rownames(res.DEA))
        group.names <- gsub("group1|group2", "", group.names)
        res.DEA$comparison <- group.names
        file <- write.table(res.DEA, paste0(filename,".txt"), sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
    }
    return(res.DEA)
}
