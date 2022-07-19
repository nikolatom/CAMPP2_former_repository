#' @title Make Heatmap
#' @description A function for making heatmaps to showcase the difference in expression values across samples.
#' @param data A dataframe containing expression values for genes within a set of samples.
#' @param group a factor derived from metadata column selected as a sample group (e.g. diagnosis), which will be used to color and split the heatmap by. Length should be equal to ncol(data)).
#' @param prefix A prefix for the output filename.
#' @param gradient A color gradient to use for the heatmap. Default is 'blue' using the viridis library.
#' @export
#' @import ComplexHeatmap
#' @import squash
#' @import viridis
#' @import grid
#' @seealso
#' @return Heatmap
#' @examples \dontrun{
#' ...
#' }

MakeHeatmap <- function(data, groups, prefix, gradient=viridis(300, option="cividis")){

    data<-as.matrix(data)
    col_ha = HeatmapAnnotation(Groups = groups,
                               simple_anno_size = unit(0.3,'cm'), annotation_name_align = TRUE)

    png(paste0(prefix,"_Heatmap.png"),width = 30, height = 20, units = "cm", res=1200)

    map <- makecmap(round(min(data)):round(max(data)))
    map$colors <- viridis((length(map$breaks)-1), option="cividis")
    lgd <- list(title='Expression value', at=1:length(map$breaks), labels=map$breaks, col_fun=map$colors)

    ht <- Heatmap(scale(data,scale=TRUE), col = gradient, heatmap_legend_param = lgd,
                  row_names_side = "left", row_names_gp=gpar(cex=0.6), column_names_gp = gpar(cex=0),
                  clustering_method_rows = "ward.D", column_title="Samples", row_title="Names",
                  top_annotation=col_ha, column_km=length(unique(groups)))

    draw(ht)
    dev.off()
}
