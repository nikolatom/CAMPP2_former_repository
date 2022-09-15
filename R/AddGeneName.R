#' @title Add HUGO Gene Name
#' @description A function for adding HUGO IDs to existing dataframe based on ENSEMBL gene IDs using biomaRt.
#' @param data a data matrix containing a column of stable gene IDs.
#' @param ensembl.version a number specifying version of ENSEMBL database is used when transforming ENSEMBL IDs into HUGO IDs using biomaRt. Default is 104.
#' @param ensembl.id.column.name a column name including ENSEMBL IDs. Default is "name" (as present in the output from DEA analysis done on limma).
#' @export
#' @import biomaRt
#' @seealso
#' @return the input dataframe merged with a column of HUGO IDs.
#' @examples \dontrun{
#' ...
#' }

AddGeneName <- function(data, ensembl.version = 104, ensembl.id.column.name="name") {

    ensembl <- useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl',version = ensembl.version)

    ensembl_id <- as.data.frame(data$ensembl.id.column.name)
    colnames(ensembl_id) <- "Gene stable ID"

    gene_names <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                        mart <- ensembl,
                        values = ensembl_id,
                        filters = 'ensembl_gene_id',
                        uniqueRows = TRUE,
                        bmHeader = T)

    names(data)[names(data) == ensembl.id.column.name] <- "Gene stable ID"
    data_merged <- merge(data, gene_names, by="Gene stable ID", all = TRUE)

    colnames(data_merged)[1] <- ensembl.id.column.name
    data_merged$Ensembl_ID <- data_merged$name
    colnames(data_merged)[which(names(data_merged) == "Gene name")] <- "HUGO_ID"

    print(paste0(" - You have chosen to use HUGO IDs for annotation. *Note that ",sum(is.na(data_merged))," genes will be lost!"))

    return(data_merged)
}
