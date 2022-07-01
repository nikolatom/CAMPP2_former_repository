#' @title Make Volcano Plot
#' @description A function for making volcano plots from DEA results.
#' @param data A data matrix from RunDEA containing genes, AveExpr, logFCs,p.value, adj.p.value ect.
#' @param prefix A prefix for the output filename.
#' @param FDR A cut-off value for  corrected p.value
#' @param logFC A cut-off value for log fold change
#' @export
#' @import ggplot2
#' @import ggrepel
#' @seealso
#' @return Volcano plot
#' @examples \dontrun{
#' ...
#' }

MakeVolcano <- function(data, prefix, logFC, FDR) {

    data <- data[order(data$P.Value),]

    data$diffexpressed <- "NO"
    data$diffexpressed[data$logFC > logFC & data$P.Value < FDR] <- "UP"
    data$diffexpressed[data$logFC < -logFC & data$P.Value < FDR] <- "DOWN"

    data$DEAlabel <- ""
    data[1:15,]$DEAlabel[data[1:15,]$diffexpressed != "NO"] <- data[1:15,]$"HUGO_ID"[data[1:15,]$diffexpressed != "NO"]

    ggplot(data, aes(x=logFC, y=-log10(P.Value), col=data$diffexpressed, label=data$DEAlabel)) +
        geom_point(shape=1) +
        theme_minimal() +
        geom_text_repel() +
        scale_color_manual(values=c("cyan3", "indianred2")) +
        geom_vline(xintercept=c(-logFC, logFC), col="black") +
        geom_hline(yintercept=-log10(FDR), col="black")

    ggsave(paste0(prefix, "_VolcanoPlot.png"), dpi = 300, width = 10, height = 8)
}
