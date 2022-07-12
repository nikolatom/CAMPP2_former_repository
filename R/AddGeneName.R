#' @title Add Gene Name function
#' @description A function for adding a column of HUGO IDs to a dataframe with a column of Ensemble IDs using biomaRt.
#' @param data a data matrix containing a column of stable gene IDs (column must be called 'name').
#' @param ensembl.version This arugment specifies which ensembl database to use when transforming ensemble IDs into HUGO IDs using biomaRt. The argument should be specified as a number.
#' @export
#' @import biomaRt
#' @seealso
#' @return the input dataframe merged with a column of HUGO IDs.
#' @examples \dontrun{
#' ...
#' }

AddGeneName <- function(data, ensembl.version) {

    ensembl <- useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl',version = ensembl.version)

    ensembl_id <- as.data.frame(data$name)
    colnames(ensembl_id) <- "Gene stable ID"

    gene_names <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                        mart <- ensembl,
                        values = ensembl_id,
                        filters = 'ensembl_gene_id',
                        uniqueRows = TRUE,
                        bmHeader = T)

    names(data)[names(data) == 'name'] <- "Gene stable ID"
    data_merged <- merge(gene_names, data, by="Gene stable ID", all = TRUE)

    colnames(data_merged)[1] <- "Ensembl_ID"
    colnames(data_merged)[2] <- "HUGO_ID"

    return(data_merged)
}
