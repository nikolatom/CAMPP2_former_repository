#' @title Download string database
#' @description Function to download string database
#' @param my.geneIDs a string specifying type of gene ID to use to match IDs in STRING database. Options are: ensembl_peptide_id, hgnc_symbol, ensembl_gene_id, ensembl_transcript_id or uniprotswissprot.
#' @param my.version Version of STRING database to use. Default is version 11.
#' @export
#' @import igraph
#' @import biomaRt
#' @import multiMiR
#' @import devtools
#' @seealso
#' @return
#' @examples \dontrun{
#' ...
#' }

DownloadPPInt <- function(my.geneIDs, my.version = "11.0") {
    approvedGeneIDs <- c("ensembl_peptide_id", "hgnc_symbol","ensembl_gene_id","ensembl_transcript_id", "uniprotswissprot")
    setwd("..")
    file <- try(load("stringDB.Rdata"))
    setwd("./")

    if (class(file) == "try-error") {
        print("\nDownloading and preparing string database for protein-protein interactions - this may take a couple of minutes!\n")
        download.file(paste0("https://stringdb-static.org/download/protein.links.v",my.version,"/9606.protein.links.v",my.version,".txt.gz"), paste0("9606.protein.links.v", my.version, ".txt.gz"))
        stringDB <- data.table::fread(paste0("9606.protein.links.v", my.version, ".txt.gz"))

        stringDB$protein1 <- gsub("9606.", "", stringDB$protein1)
        stringDB$protein2 <- gsub("9606.", "", stringDB$protein2)

        if (my.geneIDs %in% approvedGeneIDs[-1]) {
            uqprot <- unique(c(stringDB$protein1, stringDB$protein2))
            uqprot <- split(uqprot, ceiling(seq_along(uqprot)/10000))
            print("Calling IDs from BiomaRt")
            ensemblMT <- useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
            #map <- unique(getBM(attributes = c("ensembl_peptide_id", my.geneIDs), filters = "ensembl_peptide_id", values = uqprot, mart = ensemblMT))
            map <- lapply(uqprot, function(x) getBM(attributes = c("ensembl_peptide_id", my.geneIDs), filters = "ensembl_peptide_id", values = x, mart = ensemblMT))
            map <- do.call("rbind", map)
            colnames(map) <- c("protein1", "ID1")
            stringDB <- merge(stringDB, map, by = "protein1")
            colnames(map) <- c("protein2", "ID2")
            stringDB <- merge(stringDB, map, by = "protein2")
            stringDB <- stringDB[,-c(1,2)]
        }

        stringDB <- stringDB[stringDB$combined_score > as.numeric(quantile(stringDB$combined_score)[2]),]
        setwd("..")
        save(stringDB, file = "stringDB.Rdata")
        setwd("./")
    }

    return(stringDB)
}

