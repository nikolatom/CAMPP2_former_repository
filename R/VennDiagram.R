#' @title Venn Diagram
#' @description A function for making Venn diagrams to visualize the size of
#' intersections of the features between the groups
#' (e.g., differentially expressed genes for each subtype)
#' @param prefix A prefix for the output filename.
#' @param groups.feature.list A list of character vectors (each character
#' vector named by its group name, e.g., "LumA") including feature names
#' characteristic for each group of interest (e.g., subtype).
#' @export
#' @import VennDiagram
#' @import viridis
#' @seealso
#' @return Venn Diagram
#' @examples {
#' MakeVennDiagram("test_Venn", campp2_brca_1_DEA_HUGO_features_per_group)
#' }

MakeVennDiagram <- function(prefix,groups.feature.list) {

    n_col <- length(groups.feature.list)
    palette_col <- viridis_pal(alpha = 0.8,option = "viridis")(n_col)

    venn.diagram(groups.feature.list,filename = paste0(prefix,"_VennDiagram.png"),category.names = names(groups.feature.list),
                 output = TRUE,main = "Venn Diagram",main.cex = 1.5,main.fontface = "bold",
                 main.fontfamily = "sans",lwd = 3,lty = "blank",fill = palette_col,cex = 2,
                 fontface = "bold",fontfamily = "sans",cat.cex = 1.5,cat.fontface = "bold",
                 cat.fontfamily = "sans")

}
