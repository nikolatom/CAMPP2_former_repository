#' @title Make Heatmap
#' @description Function for making heatmaps.
#' @param data A dataframe with counts for differential expressed/abundant features.
#' @param gradient A color gradient to use for heatmap.
#' @param side.colors A color pallet for groups full length.
#' @param colors A color pallet for groups at levels.
#' @param group A vector of groups to color by.
#' @param filename A name of output plot.
#' @param range A vector containing the smallest and biggest logFC values from the input data frame.
#' @export
#' @import heatmap.plus
#' @import squash
#' @import viridis
#' @import ggplot2
#' @seealso
#' @return heatmaps
#' @examples \dontrun{
#' ...
#' }


MakeHeatmap <-  function(data, gradient, side.colors, colors, group, filename, range, prefix) {
    pdf(paste0(filename,"_heatmap.pdf"), height = 13, width=11)
    heatmap.plus(as.matrix(scale(data, scale = FALSE)), col=gradient, Rowv=NULL, hclustfun=function(d) hclust(d, method="ward.D2"), trace="none", labRow=rownames(data), labCol='', ColSideColors=cbind(side.colors, rep("white", length(group))), margins = c(14,8), cexCol=1, cexRow = 0.8)
    map <- makecmap(range[1]:range[2])
    map$colors <- viridis((length(map$breaks)-1), option="cividis")
    hkey(map, x = 0, y = 0, title = "LogFC", stretch = 2)
    legend("topleft", legend = as.character(levels(as.factor(group))), fill=colors, cex = 0.7)
    dev.off()
}
