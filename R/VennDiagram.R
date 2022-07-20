#' @title Venn Diagram
#' @description A function for making Venn diagrams to visualize the intersection size of significantly differentially expressed genes between 2 or more subtypes. The genes to compare is gathered from the results provided by DEA using CAMPP2.
#' @param prefix A prefix for the output filename.
#' @param sets_list A list of sets containing the gene names for each subtype to compare.
#' @param names_sets A vector of the corresponding names of each subtype in sets_list. names_sets must be same length as sets_list.
#' @export
#' @import VennDiagram
#' @import viridis
#' @seealso
#' @return Venn Diagram
#' @examples \dontrun{
#' ...
#' }

MakeVennDiagram <- function(prefix,sets_list,names_sets) {

    names_sets <- gsub("-","",gsub("healthy","",names_sets)) #As the input genes are based on DEA, the names_sets are actually comparisons, e.g. CN_high-healthy. Here we are removing everything but the subtype names to get a prettier output
    names(sets_list) <- names_sets

    n_col <- length(sets_list)
    palette_col <- viridis_pal(alpha = 0.8,option = "viridis")(n_col)

    venn.diagram(sets_list,filename = paste0(prefix,"_VennDiagram.png"),category.names = names_sets,
                 output = TRUE,main = "Venn Diagram",main.cex = 1.5,main.fontface = "bold",
                 main.fontfamily = "sans",lwd = 3,lty = "blank",fill = palette_col,cex = 2,
                 fontface = "bold",fontfamily = "sans",cat.cex = 1.5,cat.fontface = "bold",
                 cat.fontfamily = "sans")

}
