#' @title Export DEA results into .txt
#' @description A function for reformating the output from DEA and exporting the results into a .txt file.
#' @param res.DEA A list of data frames from DEA_feature_apply containing gene counts, p-values, FDRs and logFCs.
#' @param filename The name of the output file.
#' @export
#' @seealso
#' @return DEA results as .txt table
#' @examples \dontrun{
#' ...
#' }

ExportDEA <- function(res.DEA, filename) {
    if (is.null(res.DEA)) {
        cat("\nDEA yielded no results. Aren't your logFC or FDR cut-offs too strict?\n")
    } else {
        DEA.out <- do.call(rbind, unlist(res.DEA, recursive=FALSE))  ##this should be done also outside of this function
        group.names <- gsub("1[.](.*)|2[.](.*)", "", rownames(DEA.out))
        group.names <- gsub("group", "", group.names) ###removing "group" from comparison column
        DEA.out$comparison <- group.names
        write.table(DEA.out, paste0(filename,".txt"), sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        return(DEA.out)
    }
}
