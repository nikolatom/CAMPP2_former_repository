#' @title Make Volcano Plot
#' @description A function for making volcano plots from DEA results.
#' @param data a data matrix from RunDEA containing genes, AveExpr, logFCs,p.value, adj.p.value ect.
#' @export
#' @import ggplot2
#' @import ggrepel
#' @seealso
#' @return a volcano plot showcasing significantly differentially expressed genes.
#' @examples \dontrun{
#' ...
#' }

MakeVolcano <- function(data) {

    data <- data[order(data$adj.P.Val),]

    data$diffexpressed <- "NO"
    data$diffexpressed[data$logFC > 1.4 & data$adj.P.Val < 0.05] <- "UP"
    data$diffexpressed[data$logFC < -1.4 & data$adj.P.Val < 0.05] <- "DOWN"

    data$DEAlabel <- NA
    data[1:10,]$DEAlabel[data[1:10,]$diffexpressed != "NO"] <- data[1:10,]$"Gene name"[data[1:10,]$diffexpressed != "NO"]

    ggplot(data, aes(x=logFC, y=-log10(adj.P.Val), col=data$diffexpressed, label=data$DEAlabel)) +
        geom_point() +
        theme_minimal() +
        geom_text_repel() +
        scale_color_manual(values=c("blue", "black", "red")) +
        geom_vline(xintercept=c(-1.4, 1.4), col="black") +
        geom_hline(yintercept=-log10(0.05), col="black")
}
