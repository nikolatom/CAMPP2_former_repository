#' @title Make Volcano Plot of DEA features
#' @description A function for making volcano plots from DEA results.
#' The input is typically an output from RunDEA (re-formatted matrix of
#' differential expression/abundance results from limma, please, check RunDEA
#' for more details). By default, top 10 up- and down-regulated features
#' will be labeled on the Volcano plot. Important: Features are selected from
#' the table with all the groups' comparisons.
#' @param data A data frame from RunDEA containing genes and their
#' corresponding AveExpr, logFC, p.value, adj.p.value ect. Mandatory columns
#' are "logFC", "adj.P.Val" and "HUGO_ID." HUGO_ID could be obtained using
#' AddGeneName function.
#' @param prefix a prefix for the output file name.
#' @param cutoff.FDR a numeric value for adj.P.Val cutoff. Default = 0.01
#' @param cutoff.logFC a numeric value for logFC cutoff. Default = 1
#' @param n.labeled.features a number of the most up-/down- regulated features
#' labeled in the volcano plot. Default = 10
#' @export
#' @import ggplot2
#' @import ggrepel
#' @return Volcano plot of differentially expressed features saved into .png
#' @examples {
#' ##The input is made by example script from AddGeneName function
#' MakeVolcano(campp2_brca_1_DEA_HUGO, prefix = "test_volcano",
#' cutoff.logFC = 1, cutoff.FDR = 0.05, n.labeled.features = 10)
#' }

MakeVolcano <- function(data, prefix, cutoff.logFC = 1, cutoff.FDR = 0.01, n.labeled.features = 10) {

    data <- data[order(data$logFC),]  ##sort according to the logFC from the lowest value

    data$diffexpressed <- "NO"
    data$diffexpressed[data$logFC > cutoff.logFC & data$adj.P.Val < cutoff.FDR] <- "UP" ##this is based on settings specific to volcano plot (might be different to DEA settings)
    data$diffexpressed[data$logFC < -cutoff.logFC & data$adj.P.Val < cutoff.FDR] <- "DOWN" ##this is based on settings specific to volcano plot (might be different to DEA settings)

    ##Set default DEAlabel value
    data$DEAlabel <- ""

    ##Label xx most important DOWN-regulated genes
    data[1:n.labeled.features,]$DEAlabel <- head(data,n=n.labeled.features)$"HUGO_ID" ##assign labels to xx most significantly differentially expressed features

    ##Check how many of the features from top x down-regulated genes are not significant based on logFC and FDR cutoffs
    if(sum(data[1:n.labeled.features,]$diffexpressed == "DOWN") != n.labeled.features){
        num.nonsig <- n.labeled.features - (sum(data[1:n.labeled.features,]$diffexpressed == "DOWN"))
        print(paste0("Warning: ",num.nonsig, " out of ",n.labeled.features," downregulated features labeled on the Volcano plot are non-signifficant according to the FDR and logFC cutoffs."))
    } else {
        print(paste0("All ",n.labeled.features, " down-regulated features labeled on the Volcano plot are signifficantly differentially expressed."))
    }

    ##Label xx most important UP-regulated genes (taken from the end of the sorted data)
    data[(nrow(data)-(n.labeled.features-1)):(nrow(data)),]$DEAlabel <- tail(data,n=n.labeled.features)$"HUGO_ID" ##assign labels to xx most significantly differentially expressed features

    ##Check how many of the features from top x up-regulated genes are not significant based on logFC and FDR cutoffs
    if(sum(data[(nrow(data)-(n.labeled.features-1)):(nrow(data)),]$diffexpressed == "UP") != n.labeled.features){
        num.nonsig <- n.labeled.features - (sum(data[1:n.labeled.features,]$diffexpressed == "UP"))
        print(paste0("Warning: ",num.nonsig, " out of ",n.labeled.features," upregulated features labeled on the Volcano plot are non-signifficant according to the FDR and logFC cutoffs."))
    } else {
        print(paste0("All ",n.labeled.features, " up-regulated features labeled on the Volcano plot are signifficantly differentially expressed."))
    }

    ##Print volcano plot using ggplot
    ggplot(data, aes(x=logFC, y=-log10(adj.P.Val), col=diffexpressed, label=DEAlabel)) +
        geom_point(shape=1) +
        theme_minimal() +
        geom_text_repel() +
        scale_color_manual(values=c("cyan3", "indianred2","black")) +
        geom_vline(xintercept=c(-cutoff.logFC, cutoff.logFC), col="black") +
        geom_hline(yintercept=-log10(cutoff.FDR), col="black")
    ##Save the plot
    ggsave(paste0(prefix, "_VolcanoPlot.png"), dpi = 300, width = 10, height = 8)
}
