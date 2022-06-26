#' @title Heatmaps
#' @description Function for making heatmaps
#' @param my.DE a dataframe with counts for differential expressed/abundant features
#' @param my.gradient a color gradient to use for heatmap
#' @param my.colorshm a color pallet for groups full length
#' @param my.colors a color pallet for groups at levels
#' @param my.group a vector of groups to color by
#' @param my.filename a name of output plot
#' @param my.range a vector containing the smallest and biggest logFCs from differential expression/abundance analysis.
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


MakeHeatmap <-  function(my.DE, my.gradient, my.colorshm, my.colors, my.group, my.filename, my.range) {
    pdf(paste0(my.filename,"_heatmap.pdf"), height = 13, width=11)
    heatmap.plus(as.matrix(scale(my.DE, scale = FALSE)), col=my.gradient, Rowv=NULL, hclustfun=function(d) hclust(d, method="ward.D2"), trace="none", labRow=rownames(my.DE), labCol='', ColSideColors=cbind(my.colorshm, rep("white", length(my.group))), margins = c(14,8), cexCol=1, cexRow = 0.8)
    map <- makecmap(my.range[1]:my.range[2])
    map$colors <- viridis((length(map$breaks)-1), option="cividis")
    hkey(map, x = 0, y = 0, title = "LogFC", stretch = 2)
    legend("topleft", legend = as.character(levels(as.factor(my.group))), fill=my.colors, cex = 0.7)
    dev.off()
}
