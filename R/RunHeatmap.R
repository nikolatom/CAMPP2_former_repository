#' @title RunHeatmap
#' @description A function for making heatmaps to showcase the difference in
#' expression/abundance of features across samples and sample groups and to
#' visualise relations between them.
#' @param feature.counts a feature count matrix (ideally normalized and batch corrected)
#' from "seq", "array", "ms" or "other" technology (with feature IDs as row
#' names and sample IDs as columns). It's recommended to import feature counts
#' using function "import_counts". \
#' @param DEA.out a data frame from RunDEA containing genes and their
#' corresponding AveExpr, logFC, p.value, adj.p.value ect. Mandatory columns
#' are "logFC".
#' @param group a factor specifying group for each sample (e.g. could be
#' represented by a column from a metadata file). Group's information
#' will be used in clustering the samples.
#' @param heatmap.size an argument specifying how many genes will be selected
#' to be plotted on the heatmap if plot.heatmap is TRUE. The input must be
#' specified as an even number.
#' @param viridis.palette a character vector specifying color gradient to use for
#' the heatmap. Default option for viridis color palettes is 'turbo'. For more
#' options, see viridis vignette.
#' @param plot.heatmap an argument defining which data will be used for
#' the selection of the top x features to be plotted on the heatmap.
#' Options are:"DEA", "LASSO", "EN" or "Consensus".
#' @param prefix a character vector defining prefix of output file name.
#' @export
#' @import ComplexHeatmap
#' @import squash
#' @import viridisLite
#' @import grid
#' @seealso
#' @return
#' @examples \dontrun{
#' RunHeatmap(campp2_brca_1_batchCorrected, campp2_brca_1_DEA_HUGO,
#' campp2_brca_1_meta$subtype, heatmap.size=30, viridis.palette="turbo",
#' plot.heatmap="DEA", "test_heatmap")
#' }

RunHeatmap<-function(feature.counts, DEA.out, group, heatmap.size=30, viridis.palette="turbo", plot.heatmap="DEA", prefix){

    ##Selection of the features is done based on DEA output including all groups' comparisons.
    if (plot.heatmap %in% c("DEA")) {
        signif_features <- rbind(head(DEA.out[order(DEA.out$logFC),],heatmap.size/2),tail(DEA.out[order(DEA.out$logFC),],heatmap.size/2))  ##Select top x DEA features from the top/bottom for the sorted DEA output
        signif_features <- signif_features$name
        feature.counts <- feature.counts[rownames(feature.counts) %in% signif_features,]
        print(paste0("Top ", heatmap.size, " DEA features will be selected for heatplot based on their logFC values. Please, consider that some of the TOP DEA features might be present multiple times DEA output due to the multiple groups comparisons - on the heatmap, gene will be projected only once. "))
        print(paste0("Number of duplicates in the top ", heatmap.size, " DEA features is: ", (length(signif_features)-nrow(feature.counts))))
    } else if (plot.heatmap == "Consensus") {
                feature.counts <- feature.counts[rownames(feature.counts) %in% as.character(consensus$name),]
                print(paste0("\n Top ", heatmap.size, " consensual (DEA vs LASSO/EN) features will be considered in the heatplot. DEA features are selected based on their logFC values from the DEA output including all the comparisons."))
    } else if (plot.heatmap %in% c("EN", "LASSO")) {
                feature.counts <- feature.counts[rownames(feature.counts) %in% as.character(VarsSelect$LASSO.Var.Select),]
                print(paste0("\n Top ", heatmap.size, " LASSO/EN features will be considered in the heatplot."))
    } else {
                stop("You have specified argument plot.heatmap which will produce a heatmap with results from either DEA analysis, LASSO/Elastic-Net Regression or the Consensus of these. Input to argument plot.heatmap must be a string specifying which results to plot, options are: DEA, LASSO, EN, Consensus")
    }

    # Heatmap as pdf
    MakeHeatmap(feature.counts, group, prefix, viridis.palette, "Feature counts")  ##"Feature counts" might become optional

}
