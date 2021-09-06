#' @title Plot interaction network
#' @description A function plotting interaction networks
#' @param  my.trimmed.list a list of trimmed networks, i.e. output of the function "TrimWriteInt".
#' @export
#' @import multiMiR
#' @import devtools
#' @seealso
#' @return
#' @examples \dontrun{
#' ...
#' }
#'


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
