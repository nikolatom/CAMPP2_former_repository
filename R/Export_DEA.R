#' @title Reformat and Export DEA results into .txt
#' @description A function for re-formatting the output from DEA generated based
#' on the workflow: (runDEA - DEAFeatureApply - DEAFeature - limma); and
#' exporting the results into a .txt file.
#' @param res.DEA a list of data frames containing results from different group
#' comparisons in limma format.
#' @param prefix a character vector defining prefix of output file name.
#' @export
#' @return
#' 1) a data frame of DEA results
#' 2) a .txt table of exported DEA results
#' @examples {
#' campp2_brca_1_DEA_out<-ExportDEA(res.DEA = campp2_brca_1_DEA$res.DEA,
#' prefix="test")
#' }

ExportDEA <- function(res.DEA, prefix) {
    if (is.null(res.DEA)) {
        cat("\nDEA yielded no results. Please, check if your logFC or FDR cut-offs are not too strict.\n")
    } else {
        DEA.out <- do.call(rbind, unlist(res.DEA, recursive=FALSE))
        group.names <- gsub("1[.](.*)|2[.](.*)", "", rownames(DEA.out))
        group.names <- gsub("group", "", group.names) ###removing "group" from comparison column
        DEA.out$comparison <- group.names
        write.table(DEA.out, paste0(prefix,".txt"), sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        return(DEA.out)
    }
}
