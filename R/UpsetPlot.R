#' @title Upset Plot
#' @description A function for making Upset plots to visualize the size of
#' intersections of the features present each group
#' (e.g., differentially expressed genes for each subtype)
#' @param prefix A prefix for the output filename.
#' @param groups.feature.list A list of character vectors (each character
#' vector named by its group name, e.g., "LumA") including feature names
#' characteristic for each group of interest (e.g., subtype).
#' @export
#' @import ComplexHeatmap
#' @import viridis
#' @seealso
#' @return an Upset plot showcasing the size of intersections of the features
#' between the sample groups
#' @examples \dontrun{
#' MakeUpset("testUpset",campp2_brca_1_DEA_HUGO_features_per_group)
#' }

MakeUpset <- function(prefix,groups.feature.list) {

    png(file = paste0(prefix,"_UpsetPlot.png"))

    comb_mat <- make_comb_mat(groups.feature.list)

    upset_plot <- UpSet(comb_mat,set_order = names(groups.feature.list),pt_size = unit(2,"mm"),
                lwd = 1,height = unit(4,"cm"),comb_col = rainbow(length(groups.feature.list))[comb_degree(comb_mat)],
                top_annotation = upset_top_annotation(comb_mat,height = unit(12,"cm"),
                ylim = c(0,max(comb_size(comb_mat))),bar_width = 0.7,
                axis_param = list(side = "left",at = seq(from = 0,to = max(comb_size(comb_mat)),
                by = 500)),annotation_name_side = "left",annotation_name_gp = gpar(cex = 1),
                annotation_name_offset = unit(1.5,"cm")),
                right_annotation = upset_right_annotation(comb_mat,width = unit(3,"cm"),
                gp = gpar(fill = "darkseagreen"),axis_param = list(at = seq(from = 0,
                to = max(set_size(comb_mat)),by = 2000)),annotation_name_offset = unit(1,"cm")),
                row_names_gp = gpar(fontsize = 10))

    draw_upset <- draw(upset_plot)
    order <- column_order(draw_upset)
    size <- comb_size(comb_mat)

    grid.text(label = "Upset Plot",gp = gpar(fontsize = 18),y = unit(15,"cm"),
              just = "centre")

    dev.off()
}
