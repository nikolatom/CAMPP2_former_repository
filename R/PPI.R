#' @title PROTEIN-PROTEIN interactions
#' @description A Function which extracts the protein-protein interactions where each protein (gene) is differentially abundant (expressed).
#' @param my.DB = An output of the function "DownloadPPInt", e.g. a curated string database.
#' @param my.Geneset = A dataframe with results of differential expression analysis.
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

