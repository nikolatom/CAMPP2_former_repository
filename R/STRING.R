# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Function which downloads the STRING database.
# Take arguments:
# my.geneIDs = String specifying type of gene ID to use to match IDs in STRING database. Options are: ensembl_peptide_id, hgnc_symbol, ensembl_gene_id, ensembl_transcript_id or uniprotswissprot.
# my.version = Version of STRING database to use. Default is version 11.
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#' @title STRING DATABASE
#' @description FUNCTION TO OBTAIN STRING DATABASE AND ANALYZE PROTEIN-PROTEIN INTERACTIONS
#' @param my.geneIDs = String specifying type of gene ID to use to match IDs in STRING database. Options are: ensembl_peptide_id, hgnc_symbol, ensembl_gene_id, ensembl_transcript_id or uniprotswissprot.
#' @param my.version = Version of STRING database to use. Default is version 11.
#' @param my.DB = output of the function "DownloadPPInt", e.g. a curated string database.
#' @param my.Geneset = a dataframe with results of differential expression analysis.
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
    setwd("Results/")

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
        setwd("Results/")
    }

    return(stringDB)
}

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Function which subsets a dataframe of differential expression analysis into a list of dataframes by comparison.
# Take arguments:
# my.data = a dataframe with results of differential expression analysis.
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



DFtoList <- function(my.data) {

    my.data$name <- gsub("_", "-", my.data$name)
    comparisons <- levels(as.factor(my.data$comparison))
    df.list <- list()

    for (idx in 1:length(comparisons)) {
        df <- my.data[my.data$comparison == as.character(comparisons[idx]),]
        df.list[[idx]] <- df[,-c(2:4,6)]
    }

    names(df.list) <- comparisons

    df.lengths <- as.numeric(unlist(lapply(df.list, function(x) nrow(x))))
    df.lengths.min <- which(df.lengths < 2)
    if (length(df.lengths.min) > 0) {
        df.list <- df.list[-df.lengths.min]
    }
    return(df.list)
}




# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Function which extracts the protein-protein interactions where each protein (gene) is differentially abundant (expressed).
# Take arguments:
# my.DB = output of the function "DownloadPPInt", e.g. a curated string database.
# my.Geneset = a dataframe with results of differential expression analysis.
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


GetDeaPPInt <- function(my.DB, my.Genedat) {

    df.list <- DFtoList(my.Genedat)
    nodes.list <- list()

    for (idx in 1:length(df.list)) {
        # Merging StringDB and datasets
        df <- df.list[[idx]]
        df$ID1 <- df$name
        nodes <- merge(my.DB, df, by = "ID1")
        df$ID2 <- df$name
        nodes <- unique(merge(nodes, df, by = "ID2"))
        nodes <- nodes[order(nodes$combined_score, decreasing = TRUE),]
        df <- as.data.frame(nodes[, c(2,1,3,4,5,7,9,10,12,13)])
        colnames(df) <- c("node1", "node2", "score", "logFC.node1", "fdr.node1", "dir.node1", "logFC.node2", "fdr.node2", "dir.node2", "comparison")
        df <- df[!duplicated(apply(df,1,function(x) paste(sort(x),collapse=''))),]
        nodes.list[[idx]] <- df
    }

    names(nodes.list) <- names(df.list)
    return(nodes.list)
}

