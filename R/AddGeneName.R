#' @title Add Gene Name function
#' @description A function for transforming stable gene IDs into HUGO IDs using biomaRt.
#' @param data a data matrix containing a column of stable gene IDs (column must be called 'name').
#' @export
#' @import biomaRt
#' @seealso
#' @return the input data matrix merged with a gene name column.
#' @examples \dontrun{
#' ...
#' }

AddGeneName <- function(data) {

    ensembl104 <- useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl',version = 104)

    ensembl_id <- as.data.frame(data$name)
    colnames(ensembl_id) <- "Gene stable ID"

    gene_names <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                        mart = ensembl104,
                        values = ensembl_id,
                        filters = 'ensembl_gene_id',
                        uniqueRows = TRUE,
                        bmHeader = T)

    names(data)[names(data) == 'name'] <- "Gene stable ID"
    data_merged <- merge(gene_names, data, by="Gene stable ID", all = TRUE)

    colnames(data_merged)[1] <- "Ensembl ID"
    colnames(data_merged)[2] <- "HUGO ID"

    return(data_merged)
}
