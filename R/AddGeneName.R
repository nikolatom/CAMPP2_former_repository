#' @title Add HUGO Gene Name
#' @description A function for adding HUGO IDs ("external_gene_name") to
#' existing data frame based on ENSEMBL gene IDs using biomaRt. The input is
#' typically an output from RunDEA (re-formatted matrix of
#' differential expression/abundance results from limma, please, check RunDEA
#' for more details). The output is a data frame with a column containing HUGO IDs
#' @param data a data frame containing a column of ENSEMBL IDs.
#' @param ensembl.version a number specifying version of ENSEMBL database used
#' for obtaining HUGO IDs using biomaRt. Default is 104.
#' @param ensembl.id.column.name a column name including ENSEMBL IDs. Default
#' is "name" (based on limma output format).
#' @export
#' @import biomaRt
#' @return a data frame with a column of HUGO IDs.
#' @examples {
#' campp2_brca_1_DEA_HUGO <- AddGeneName(campp2_brca_1_DEA$DEA.out, 104)
#' }

AddGeneName <- function(data, ensembl.version = 104, ensembl.id.column.name="name") {

    ensembl <- useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl',version = ensembl.version)

    ensembl_id <- as.data.frame(data[ensembl.id.column.name])
    colnames(ensembl_id) <- "Gene stable ID"

    gene_names <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                        mart <- ensembl,
                        values = ensembl_id,
                        filters = 'ensembl_gene_id',
                        uniqueRows = TRUE,
                        bmHeader = TRUE)

    names(gene_names)[1]<-ensembl.id.column.name
    names(gene_names)[2]<-"HUGO_ID"

    data.HUGO <- merge(data, gene_names, by=ensembl.id.column.name, all = TRUE)

    data.HUGO <- data.HUGO[,c(2,3,4,5,6,1,7,8,9)]  ##reorder the columns to the original format with HUGO as last

    print(paste0(" - You have chosen to use HUGO IDs for annotation. *Note that ",sum(is.na(data.HUGO))," genes will be lost due to the lack of HUGO for a given ENSEMBL ID!"))

    return(data.HUGO)
}
