# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Function which extracts the miRNA-gene interactions where miRNAs (and potentially genes, if this data is available,) are differentially expressed.
# Take arguments:
# arg.miRset = string specifying what miRNA database to use, options are: mirtarbase (validated), targetscan (predicted) or tarscanbase (validated and predicted).
# my.miRdat = a dataframe with results of differential expression analysis (miRNAs).
# my.Geneset = a dataframe with results of differential expression analysis (genes)- by default this argument is set to NULL.
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#' @title miRNA INTERACTIONS
#' @description FUNCTION TO OBTAIN miRNA interactions
#' @param arg.miRset = string specifying what miRNA database to use, options are: mirtarbase (validated), targetscan (predicted) or tarscanbase (validated and predicted).
#' @param  my.miRdat = a dataframe with results of differential expression analysis (miRNAs).
#' @param  my.Geneset = a dataframe with results of differential expression analysis (genes)- by default this argument is set to NULL.
#' @param nodes.list1 = List of networks, i.e. output of the function "GetGeneMiRNAInt".
#' @param nodes.list2 = List of networks, i.e. output of the function "GetDeaPPInt".
#' @param my.nodes.list = List of networks, i.e. output of the function "GetGeneMiRNAInt" or ourput of "GetDeaPPInt" or output of "MatchPPGmiRInt".
#' @param my.trimmed.list = List of trimmed networks, i.e. output of the function "TrimWriteInt".
#' @export
#' @import igraph
#' @import biomaRt
#' @import multiMiR
#' @import devtools
#' @seealso
#' @return DA/DE features
#' @examples \dontrun{
#' ...
#' }
#'
GetGeneMiRNAInt <- function(arg.miRset, my.miRdat, my.Genedat = NULL) {

    miR.list <- DFtoList(my.miRdat)
    nodes.list <- list()

    for (idx in 1:length(miR.list)) {

        df <- miR.list[[idx]]
        df$node1 <- df$name


        if (arg.miRset == "mirtarbase") {
            mirs <- get_multimir(mirna = df$node1, table = "mirtarbase")@data
            if(length(mirs) > 0) {
                mirs <- unique(mirs[,3:4])
                mirs$score <- rep(1, nrow(mirs))
                colnames(mirs) <- c("node1", "node2", "score")
            } else {
                stop("Either there are no differentially expressed miRNAs or non of miRNAs the have any predicted gene targets. Try re-running the pipeline with a lower cutoff for significance (-f) or change argument -i to either targetscan or tarscanbase")
            }
        } else if (arg.miRset == "targetscan") {
            mirs <- get_multimir(mirna = df$node1, table = "targetscan")@data
            if(length(mirs) > 0) {
                mirs <- mirs[,c(3,4,7)]
                mirs$score <- abs(as.numeric(mirs$score))
                colnames(mirs) <- c("node1", "node2", "score")
            } else {
                stop("Either there are no differentially expressed miRNAs or non of miRNAs the have any predicted gene targets. Try re-running the pipeline with a lower cutoff for significance (-f) or change argument -i to either mirtarbase or tarscanbase")
            }
        } else if (arg.miRset == "tarscanbase") {
            mirsV <- get_multimir(mirna = df$node1, table = "mirtarbase")@data
            if(length(mirsV) > 0) {
                mirsV <- mirsV[, 3:4]
                mirsV$score <- rep(1, nrow(mirsV))
            }
            mirsP <- get_multimir(mirna = df$node1, table = "targetscan")@data
            if(length(mirsP) > 0) {
                mirsP <- mirsP[,c(3,4,7)]
            }
            if (length(mirsV) > 0 | length(mirsP) > 0) {
                mirs <- unique(rbind(mirsV, mirsP))
                mirs$score <- abs(as.numeric(mirs$score))

                mirs$IDs <- paste0(mirs[,1], "|", mirs[,2])
                mirs <- data.table(mirs[,-c(1,2)])
                mirs <-  data.frame(mirs[, lapply(.SD, mean), by=IDs])
                mirs <- data.frame(do.call(rbind, strsplit(mirs$IDs, "[|]")), as.numeric(mirs[,2]))
                colnames(mirs) <- c("node1", "node2", "score")
            } else {
                stop("Either there are no differentially expressed miRNAs or non of miRNAs the have any predicted gene targets. No network can be generated!")
            }
        } else {
            stop("If the argument -i is specified, the second part of this argument it must be set to either mirtarbase, targetscan or tarscanbase!")
        }

        mirs <- mirs[!duplicated(apply(mirs,1,function(x) paste(sort(x),collapse=''))),]
        mirs <- mirs[mirs$score > as.numeric(quantile(mirs$score)[3]),]
        dfmiR <- merge(df, mirs, by = "node1")
        dfmiR <- dfmiR[, c(1,7,8,2,3,5,6)]
        colnames(dfmiR) <- c("node1", "node2", "score", "logFC.node1", "fdr.node1", "dir.node1", "comparison")
        dfmiR$score <- dfmiR$score * 1000
        dfmiR <- dfmiR[-grep("^$", dfmiR$node2),]

        if(is.null(my.Genedat)) {
            dfmiR$logFC.node2 <- as.numeric(rep(0, nrow(dfmiR)))
            dfmiR$fdr.node2 <- rep(NA, nrow(dfmiR))
            dfmiR$dir.node2 <- rep(NA, nrow(dfmiR))
        }
        nodes.list[[idx]] <- dfmiR
    }

    names(nodes.list) <- names(miR.list)

    if(!is.null(my.Genedat)) {

        miRgene.list  <- list()

        gene.list <- DFtoList(my.Genedat)

        overlap <- intersect(names(gene.list), names(nodes.list))

        nodes.list <- nodes.list[names(nodes.list) %in% overlap]
        gene.list <- gene.list[names(gene.list) %in% overlap]

        for (idx in 1:length(gene.list)) {

            df <- gene.list[[idx]]
            df$node2 <- df$name

            dfmiRgene <- merge(df, nodes.list[[idx]], by = "node2")
            dfmiRgene <- dfmiRgene[, c(7,1,8,9,10,11,2,3,5,6)]
            colnames(dfmiRgene) <-c("node1", "node2", "score", "logFC.node1", "fdr.node1", "dir.node1", "logFC.node2", "fdr.node2", "dir.node2", "comparison")
            dfmiRgene <- dfmiRgene[dfmiRgene$dir.node1 != dfmiRgene$dir.node2, ]

            miRgene.list[[idx]] <- dfmiRgene
        }
        nodes.list <- miRgene.list
        names(nodes.list) <- overlap
    }
    return(nodes.list)
}


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Function for combining miRNA-gene interactions and protein-protein (gene-gene) interactions.
# Take arguments:
# nodes.list1 = List of networks, i.e. output of the function "GetGeneMiRNAInt".
# nodes.list2 = List of networks, i.e. output of the function "GetDeaPPInt".
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



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








# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Function to write out interaction lists and trim interactions for plotting
# Take arguments:
# my.nodes.list = List of networks, i.e. output of the function "GetGeneMiRNAInt" or ourput of "GetDeaPPInt" or output of "MatchPPGmiRInt".
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


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




# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Function plotting interaction networks.
# Take arguments:
# my.trimmed.list = List of trimmed networks, i.e. output of the function "TrimWriteInt".
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



PlotInt <- function(my.trimmed.list) {

    trim.lengths <- as.numeric(unlist(lapply(my.trimmed.list, function(x) nrow(x))))
    trim.lengths.min <- which(trim.lengths < 10)

    if (length(trim.lengths.min) > 0) {
        my.trimmed.list <- my.trimmed.list[-trim.lengths.min]
    }

    trimmed.names <- names(my.trimmed.list)

    for (idx in 1:length(my.trimmed.list)) {

        p.nodes <- my.trimmed.list[[idx]]$p.nodes
        p.info <- my.trimmed.list[[idx]]$p.info

        final.nodes <- as.matrix(p.nodes[, 1:2])
        degrees <- as.numeric(p.info$Freq)
        names(degrees) <- p.info$name
        meta <- data.frame(p.info$group, degrees, p.info$name, ind=1:nrow(p.info))
        my.order <- as.integer(meta[order(meta$degrees, decreasing = TRUE),]$ind)

        tiff(paste0(trimmed.names[[idx]],"_TopInteracions.tiff"), width = 14, height = 10, units = 'in', res = 600)

        #cex.nodes = 2^(log2(abs(p.info$logFC))/10)

        arcplot(final.nodes, ordering=my.order, labels=p.info$name, cex.labels=0.6,
                show.nodes=TRUE, col.nodes=p.info$calfb, bg.nodes=p.info$calfb,
                cex.nodes = (log(degrees)/2)+0.2, pch.nodes=21,
                lwd.nodes = 2, line=-0.5,
                col.arcs = hsv(0, 0, 0.2, 0.25), lwd.arcs = log2(p.nodes$score/100)*1.5)
        dev.off()
    }
}






