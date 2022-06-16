#' @title Write and trim interactions
#' @description A function to write out interaction lists and trim interactions for plotting
#' @param  my.nodes.list A list of networks, i.e. output of the function "GetGeneMiRNAInt" or ourput of "GetDeaPPInt" or output of "MatchPPGmiRInt".
#' @export
#' @import multiMiR
#' @import devtools
#' @seealso
#' @return
#' @examples \dontrun{
#' ...
#' }
#'

TrimWriteInt <- function(my.nodes.list) {

    node.lengths <- as.numeric(unlist(lapply(my.nodes.list, function(x) nrow(x))))
    node.lengths.min <- which(node.lengths < 2)

    if (length(node.lengths.min) > 0) {
        my.nodes.list <- my.nodes.list[-node.lengths.min]
    }

    set.names <- names(my.nodes.list)
    trimlist <- list()

    for (idx in 1:length(my.nodes.list)) {
        p.nodes <- my.nodes.list[[idx]]
        p.nodes$myrank <- rowMeans(apply(cbind(abs(p.nodes$logFC.node1), abs(p.nodes$logFC.node2)), 2, rank))
        p.nodes <- p.nodes[order(p.nodes$myrank, decreasing = TRUE),]
        write.table(p.nodes, paste0(set.names[[idx]],"_AllInteractions.txt"), sep = "\t", row.names=FALSE, col.names=TRUE, quote = FALSE)

        if (nrow(p.nodes) > 100) {
            miRidx <- grep("miR|let", p.nodes$node1)
            if (length(miRidx) > 0) {
                miR.nodes <- p.nodes[miRidx,]

                if (sum(miR.nodes$logFC.node2) == 0) {
                    miRidx <- which(miR.nodes$score > as.numeric(quantile(miR.nodes$score)[4]))
                    miR.nodes <- miR.nodes[miRidx,]
                }

                p.nodes <- p.nodes[-miRidx,]

                if(nrow(miR.nodes) >= 100) {
                    p.nodes <- miR.nodes[1:100,]
                } else {
                    p.nodes <- p.nodes[1:(100-nrow(miR.nodes)),]
                    p.nodes <- rbind(p.nodes, miR.nodes)
                }
            } else {
                p.nodes <- p.nodes[1:100,]
            }
        }

        p.info <- data.table(c(as.character(p.nodes$node1), as.character(p.nodes$node2)), c(p.nodes$logFC.node1, p.nodes$logFC.node2))
        colnames(p.info) <- c("name", "logFC")
        p.info<- data.frame(p.info[, .(Freq = .N), by = .(name, logFC)])
        p.info <- p.info[order(p.info$Freq, decreasing = TRUE),]
        p.info$group <- 1:nrow(p.info)
        vircls <- viridis(2, direction = -1 , end = 0.9, option = "cividis")
        #vircls <- viridis(2, end = 0.6, direction = -1, option = "magma")
        p.info$calfb <- ifelse(p.info$logFC > 0, vircls[1], "grey50")
        p.info$calfb <- ifelse(p.info$logFC < 0, vircls[2], as.character(p.info$calfb))

        myorder <- unique(as.vector(t(data.frame(p.nodes[,1:2]))))
        p.info <- p.info[match(myorder, p.info$name),]

        trimmed <- list(p.nodes, p.info)
        names(trimmed) <- c("p.nodes", "p.info")
        trimlist[[idx]] <- trimmed
    }
    names(trimlist) <- set.names
    return(trimlist)
}

