#' @title Make Volcano Plot
#' @description A function for making volcano plots from DEA results.
#' @param data A data matrix from RunDEA containing genes, AveExpr, logFCs,p.value, adj.p.value ect.
#' @export
#' @import ggplot2
#' @import ggrepel
#' @seealso
#' @return Volcano plot
#' @examples \dontrun{
#' ...
#' }

MakeVolcano <- function(data, prefix) {

    data <- data[order(data$P.Value),]

    data$diffexpressed <- "NO"
    data$diffexpressed[data$logFC > 1.0 & data$P.Value < 0.05] <- "UP"
    data$diffexpressed[data$logFC < -1.0 & data$P.Value < 0.05] <- "DOWN"

    data$DEAlabel <- ""
    data[1:10,]$DEAlabel[data[1:10,]$diffexpressed != "NO"] <- data[1:10,]$"Gene name"[data[1:10,]$diffexpressed != "NO"]

    ggplot(data, aes(x=logFC, y=-log10(P.Value), col=data$diffexpressed, label=data$DEAlabel)) +
        labs(col="Gene regulation:") +
        geom_point(shape=1) +
        theme_minimal() +
        geom_text_repel() +
        scale_color_manual(values=c("cyan3", "black", "indianred2")) +
        geom_vline(xintercept=c(-1.0, 1.0), col="black") +
        geom_hline(yintercept=-log10(0.05), col="black")
    ggsave(paste0(prefix, "_VolcanoPlot.pdf"), dpi = 300, width = 10, height = 8)
}
