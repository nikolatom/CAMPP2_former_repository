#' @title Get heatmaps colors
#' @description Function for getting heatmap colors
#' @param my.truestatus a vetcor of groups/labels (a character vector, length of ncol in the matrix to be plotted)
#' @param my.colors a vector with colors to use (a character vector with the length of the number of groups/levels).
#' @export
#' @import heatmap.plus
#' @import squash
#' @importFrom viridisLite viridis
# @import ggplot2
#' @seealso
#' @return heatmap colors
#' @examples \dontrun{
#' ...
#' }


GetColors <- function(my.truestatus, my.colors) {
    hm_col <- data.frame(status = levels(as.factor(my.truestatus)), mycolor = my.colors)
    true_status <- data.frame(status = my.truestatus)
    col <- join(true_status, hm_col)
    col$mycolor <- ifelse(is.na(col$mycolor), "black", as.character(col$mycolor))
    return(as.matrix(col$mycolor))
}
