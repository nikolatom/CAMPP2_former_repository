#' @title Upset Plot
#' @description A function for making Upset plots to visualize the intersection between 2 or more sets.
#' @param prefix A prefix for the output filename.
#' @param sets_list A list of sets to compare.
#' @param names_sets A vector of the corresponding names to each set in sets_list.
#' @param title_plot
#' @export
#' @import ComplexHeatmap
#' @import viridis
#' @seealso
#' @return Upset Plot
#' @examples \dontrun{
#' ...
#' }

MakeUpset <- function(prefix,sets_list,names_sets) {

  pdf(file = paste0(prefix,"_UpsetPlot.pdf"), width = 18,height = 8)

  names(sets_list) <- names_sets
  comb_mat_sets <- make_comb_mat(sets_list)
  n_col <- max(comb_degree(comb_mat_sets))
  palette_col <- viridis_pal(option = "viridis")(n_col)

  upset_p <- UpSet(comb_mat_sets,set_order = names(sets_list),pt_size = unit(5,"mm"),
                   lwd = 3,height = unit(4,"cm"),comb_col = palette_col[comb_degree(comb_mat_sets)],
                   top_annotation = upset_top_annotation(comb_mat_sets,height = unit(12,"cm"),
                   ylim = c(0,max(comb_size(comb_mat_sets))),bar_width = 0.7,
                   axis_param = list(side = "left",at = seq(from = 0,to = max(comb_size(comb_mat_sets)),
                   by = 500)),annotation_name_side = "left",annotation_name_gp = gpar(cex = 1),
                   annotation_name_offset = unit(1.5,"cm")),
                   right_annotation = upset_right_annotation(comb_mat_sets,width = unit(3,"cm"),
                   gp = gpar(fill = "darkseagreen"),axis_param = list(at = seq(from = 0,
                   to = max(set_size(comb_mat_sets)),by = 2000)),annotation_name_offset = unit(1.5,"cm")),
                   row_names_gp = gpar(fontsize = 18))

  draw_upset <- draw(upset_p)
  col_ord <- column_order(draw_upset)
  c_s <- comb_size(comb_mat_sets)
  decorate_annotation("intersection_size", {grid.text(c_s[col_ord],x = seq(c_s),
                                                      y = unit(c_s[col_ord],"native") + unit(2, "pt"),gp = gpar(fontsize = 12,
                                                                                                                fontface = "bold"),just = "bottom",default.units = "native")
  })

  grid.text(label = "Upset Plot",x = unit(20,"cm"),y = unit(18,"cm"),gp = gpar(fontsize = 18),just = "centre")

  dev.off()
}
