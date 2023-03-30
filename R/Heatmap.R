#' @title Make Heatmap
#' @description A function for making heatmaps to showcase the difference in
#' expression/abundance of features across samples and sample groups and to
#' visualise relations between them.
#' @param data a feature count matrix (ideally normalized and batch corrected)
#' from "seq", "array", "ms" or "other" technology (with feature IDs as row
#' names and sample IDs as columns). It's recommended to import feature counts
#' using function "import_counts".
#' @param groups a factor specifying group for each sample (e.g. could be
#' represented by a column from a metadata file). Group's information
#' will be used in clustering the samples.
#' @param prefix a character vector defining prefix of output file name.
#' @param viridis.palette a character vector specifying color gradient to use for
#' the heatmap. Default option for viridis color palettes is 'turbo'. For more
#' options, see viridis vignette.
#' @param data.type a character vector (used in the heatmap legend) describing
#' type of the data being visualized. Default is "expression/abundance".
#' @export
#' @import ComplexHeatmap
#' @import squash
#' @import viridisLite
#' @import grid
#' @return heatmap
#' @examples {
#' MakeHeatmap(data=campp2_brca_1_batchCorrected,
#' groups=campp2_brca_1_meta$subtype, prefix="test", viridis.palette = "turbo",
#' data.type = "Feature counts")
#' }

MakeHeatmap <- function(data, groups, prefix, viridis.palette = "turbo", data.type = "expression/abundance"){

    col_ha = HeatmapAnnotation(Groups = groups,
                               simple_anno_size = unit(0.3,'cm'), annotation_name_align = TRUE)

    png(paste0(prefix,"_heatmap.png"), width = 30, height = 20, units = "cm", res=1200)

    map <- makecmap(round(min(data)):round(max(data)))
    map$colors <- viridis((length(map$breaks)-1), option=viridis.palette)  ##change colors according to viris pallete
    legend <- list(title=data.type, at=1:length(map$breaks), labels=map$breaks, col_fun=map$colors)

    htmap <- Heatmap(scale(data,scale=TRUE), col = viridis(300, option=viridis.palette), heatmap_legend_param = legend,
                     row_names_side = "left", row_names_gp=gpar(cex=0.6), column_names_gp = gpar(cex=0),
                     clustering_method_rows = "ward.D", column_title="Samples", row_title="Features",
                     top_annotation=col_ha, column_km=length(unique(groups)))

    draw(htmap)
    dev.off()
}
