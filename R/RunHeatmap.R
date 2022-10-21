#' @title RunHeatmap
#' @description A function for making heatmaps to showcase the difference in
#' expression/abundance of features across samples and sample groups and to
#' visualise relations between them.
#' @param feature.counts a feature count matrix (ideally normalized and batch corrected)
#' from "seq", "array", "ms" or "other" technology (with feature IDs as row
#' names and sample IDs as columns). It's recommended to import feature counts
#' using function "import_counts".
#' @param DEA.out a data frame from RunDEA containing genes and their
#' corresponding AvgExpr, logFC, p.value, adj.p.value ect. Mandatory columns
#' are "logFC".
#' @param group a factor specifying group for each sample (e.g. could be
#' represented by a column from a metadata file). Group's information
#' will be used in clustering the samples.
#' @param viridis.palette a character vector specifying color gradient to use for
#' the heatmap. Default option for viridis color palettes is 'turbo'. For more
#' options, see viridis vignette.
#' @param plot.heatmap an argument defining which data will be used for
#' the selection of the top x features to be plotted on the heatmap.
#' Options are:"DEA", "LASSO", "EN", "Ridge" or "Consensus".
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
#' campp2_brca_1_meta$subtype, viridis.palette="turbo",
#' plot.heatmap="DEA", "test_heatmap")
#' }

RunHeatmap<-function(feature.counts, DEA.out, group, viridis.palette="turbo", plot.heatmap="DEA", prefix){

    ##Selection of the features is done based on DEA output including all groups' comparisons.
    if (plot.heatmap %in% c("DEA","Consensus","EN", "LASSO", "Ridge")) {
        feature.counts <- feature.counts[rownames(feature.counts) %in% DEA.out,]
    } else {
        stop("You have specified argument plot.heatmap which will produce a heatmap with results from either DEA analysis, LASSO/Elastic-Net Regression or the Consensus of these. Input to argument plot.heatmap must be a string specifying which results to plot, options are: DEA, LASSO, EN, Ridge, Consensus")
    }

    # Heatmap as pdf
    MakeHeatmap(feature.counts, group, prefix, viridis.palette, "Feature counts")  ##"Feature counts" might become optional

}
