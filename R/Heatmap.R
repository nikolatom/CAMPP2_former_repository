#' @title Make Heatmap
#' @description A function for making heatmaps from expression data.
#' @param data A dataframe with counts for differential expressed/abundant features.
#' @param gradient A color gradient to use for the heatmap.
#' @param groups A vector of groupss to color by.
#' @param prefix A prefix for the output filename.
#' @param range A vector containing the smallest and biggest logFC values from the input data frame.
#' @export
#' @import ComplexHeatmap
#' @import squash
#' @import viridis
#' @seealso
#' @return Heatmap
#' @examples \dontrun{
#' ...
#' }

MakeHeatmap <- function(data, gradient, groups, prefix, range){
    data<-as.matrix(data)

    col_ha = HeatmapAnnotation(Groups = group1, col = list(foo = col_fun) ,
                               simple_anno_size = unit(0.3,'cm'), annotation_name_align = TRUE)

    png(paste0(prefix,"_Heatmap.png"),width = 30, height = 20, units = "cm", res=1200)

    map <- makecmap(range[1]:range[2])
    map$colors <- viridis((length(map$breaks)-1), option="cividis")
    lgd <- list(title='LogFC', at=1:length(map$breaks), labels=map$breaks, col_fun=map$colors)

    ht <- Heatmap(scale(data,scale=TRUE), col = gradient, heatmap_legend_param = lgd,
                  row_names_side = "left", row_names_gp=gpar(cex=0.6), column_names_gp = gpar(cex=0),
                  clustering_method_rows = "ward.D", column_title="Samples", row_title="Genes",
                  top_annotation=col_ha, column_km=length(unique(groups)))

    draw(ht)
    dev.off()
}
