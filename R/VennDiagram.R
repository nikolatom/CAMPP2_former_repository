#' @title Venn Diagram
#' @description A function for making Venn Diagrams to visualize the intersection between 2 or more sets.
#' @param prefix A prefix for the output filename.
#' @param sets_list A list of sets to compare.
#' @param names_sets A vector of the corresponding names to each set in sets_list.
#' @export
#' @import VennDiagram
#' @import viridis
#' @seealso
#' @return Venn Diagram
#' @examples \dontrun{
#' ...
#' }

MakeVennDiagram <- function(sets_list,names_sets,prefix) {

    # Assign names to list of sets
    names(sets_list) <- names_sets

    # Define colors to use in Venn diagram
    n_col <- length(sets_list)
    palette_col <- viridis_pal(alpha = 0.8,option = "viridis")(n_col)

    #Create Venn diagram and save it as pdf in plots directory of TCGA project
    venn.diagram(sets_list,filename = paste0(prefix,"_VennDiagram.png"),category.names = names_sets,
                 output = TRUE,main = "Venn Diagram",main.cex = 1.5,main.fontface = "bold",
                 main.fontfamily = "sans",lwd = 3,lty = "blank",fill = palette_col,cex = 2,
                 fontface = "bold",fontfamily = "sans",cat.cex = 1.5,cat.fontface = "bold",
                 cat.fontfamily = "sans")

}
