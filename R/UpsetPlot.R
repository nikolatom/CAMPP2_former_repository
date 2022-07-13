#' @title Upset Plot
#' @description A function for making Upset plots to visualize the intersection size of significantly differentially expressed genes between 2 or more subtypes. The genes to compare is gathered from the results provided by DEA using CAMPP2.
#' @param prefix A prefix for the output filename.
#' @param sets_list A list of sets containing the gene names for each subtype to compare.
#' @param names_sets A vector of the corresponding names of each subtype in sets_list. names_sets must be same length as sets_list.
#' @export
#' @import ComplexHeatmap
#' @import viridis
#' @seealso
#' @return Upset plot
#' @examples \dontrun{
#' ...
#' }

MakeUpset <- function(prefix,sets_list,names_sets) {

    png(file = paste0(prefix,"_UpsetPlot.png"))

    names_sets <- gsub("-","",gsub("healthy","",names_sets))
    names(sets_list) <- names_sets
    comb_mat <- make_comb_mat(sets_list)
    n_col <- max(comb_degree(comb_mat))

    upset_p <- UpSet(comb_mat,set_order = names(sets_list),pt_size = unit(2,"mm"),
                lwd = 1,height = unit(4,"cm"),comb_col = rainbow(length(sets_list))[comb_degree(comb_mat)],
                top_annotation = upset_top_annotation(comb_mat,height = unit(12,"cm"),
                ylim = c(0,max(comb_size(comb_mat))),bar_width = 0.7,
                axis_param = list(side = "left",at = seq(from = 0,to = max(comb_size(comb_mat)),
                by = 500)),annotation_name_side = "left",annotation_name_gp = gpar(cex = 1),
                annotation_name_offset = unit(1.5,"cm")),
                right_annotation = upset_right_annotation(comb_mat,width = unit(3,"cm"),
                gp = gpar(fill = "darkseagreen"),axis_param = list(at = seq(from = 0,
                to = max(set_size(comb_mat)),by = 2000)),annotation_name_offset = unit(1,"cm")),
                row_names_gp = gpar(fontsize = 10))

    draw_upset <- draw(upset_p)
    order <- column_order(draw_upset)
    size <- comb_size(comb_mat)

    grid.text(label = "Upset Plot",gp = gpar(fontsize = 18),y = unit(15,"cm"),
              just = "centre")

    dev.off()
}
