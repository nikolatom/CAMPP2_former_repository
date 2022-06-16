#' @title GetGeneMiRNAInt
#' @description A function which extracts the miRNA-gene interactions where miRNAs (and potentially genes, if this data is available,) are differentially expressed.
#' @param arg.miRset a string specifying what miRNA database to use, options are: mirtarbase (validated), targetscan (predicted) or tarscanbase (validated and predicted).
#' @param  my.miRdat a dataframe with results of differential expression analysis (miRNAs).
#' @param  my.Genedat ?
#' @export
#' @import multiMiR
#' @import devtools
#' @seealso
#' @return
#' @examples \dontrun{
#' ...
#' }


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

