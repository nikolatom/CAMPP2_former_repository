#' @title Match PP GmiR Intactions
#' @description A function for combining miRNA-gene interactions and protein-protein (gene-gene) interactions.
#' @param  nodes.list1 a list of networks, i.e. output of the function "GetGeneMiRNAInt"
#' @param  nodes.list2 a list of networks, i.e. output of the function "GetDeaPPInt"
#' @export
#' @import multiMiR
#' @import devtools
#' @seealso
#' @return
#' @examples \dontrun{
#' ...
#' }



MatchPPGmiRInt <- function(nodes.list1, nodes.list2) {

    nodes.list <- list()

    overlap <- intersect(names(nodes.list1), names(nodes.list2))

    nodes.list1 <- nodes.list1[names(nodes.list1) %in% overlap]
    nodes.list2 <- nodes.list2[names(nodes.list2) %in% overlap]

    for (idx in 1:length(nodes.list1)) {

        df1 <- nodes.list1[[idx]]
        df2 <- nodes.list2[[idx]]

        nodes <- rbind(df1[df1$dir.node1 == "up",], df1[df1$dir.node1 == "down",], df2[df2$dir.node1 == "up",], df2[df2$dir.node1 == "down",])
        nodes.list[[idx]] <- nodes
    }

    names(nodes.list) <- overlap
    return(nodes.list)
}

