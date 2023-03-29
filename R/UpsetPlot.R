#' @title Upset Plot
#' @description A function for making Upset plots to visualize the size of
#' intersections of the features present each group
#' (e.g., differentially expressed genes for each subtype).
#' @param prefix a prefix for the output filename.
#' @param groups.feature.list a list of character vectors (each character
#' vector named by its group name, e.g., "LumA") including feature names
#' characteristic for each group of interest (e.g., subtype)
#' @param label a character vector used as a label for the upset plot. Default
#' is "".
#' @param label.vertical.position a numeric value describing vertical position of
#' the label. Default is 15.0 (cm)
#' @param y.axis.by an integer describing interval used for labeling Y-axis.
#' Defaults is 100.
#' @param set.size.by an integer describing interval used for labeling the
#' set size axis. Defaults is 200.
#' @export
#' @import ComplexHeatmap
#' @import viridis
#' @return an Upset plot showcasing the size of intersections of the features
#' between the sample groups
#' @examples {
#' control.group="healthy"
#' groups.full.list <- split.data.frame(campp2_brca_1_DEA_HUGO, campp2_brca_1_DEA_HUGO$comparison)
#' campp2_brca_1_DEA_HUGO_features_per_group=list()
#' for (features in groups.full.list){
#'     campp2_brca_1_DEA_HUGO_features_per_group <- append(campp2_brca_1_DEA_HUGO_features_per_group,list(features$name))
#' }
#' names(campp2_brca_1_DEA_HUGO_features_per_group) <- gsub("-","",gsub(control.group,"",names(groups.full.list)))

#' MakeUpset("testUpSet",campp2_brca_1_DEA_HUGO_features_per_group,
#' label="test_UpSet",label.vertical.position=16.5,
#' y.axis.by=300,set.size.by=250)
#' }

MakeUpset <- function(prefix,groups.feature.list,label="",label.vertical.position=15.0,y.axis.by=100,set.size.by=200) {

    pdf(file = paste0(prefix,"_UpsetPlot.pdf"))

    comb_mat <- make_comb_mat(groups.feature.list)

    upset_plot <- UpSet(comb_mat,set_order = names(groups.feature.list),
                pt_size = unit(2,"mm"),
                lwd = 1,
                height = unit(4,"cm"),
                comb_col = rainbow(length(groups.feature.list))[comb_degree(comb_mat)],
                top_annotation = upset_top_annotation(comb_mat,height = unit(12,"cm"),
                ylim = c(0,max(comb_size(comb_mat))),
                bar_width = 0.7,
                axis_param = list(side = "left",at = seq(from = 0,to = max(comb_size(comb_mat)),
                by = y.axis.by)),
                annotation_name_side = "left",
                annotation_name_gp = gpar(cex = 1),
                annotation_name_offset = unit(1.5,"cm")),
                right_annotation = upset_right_annotation(comb_mat,width = unit(3,"cm"),
                gp = gpar(fill = "darkseagreen"),
                axis_param = list(at = seq(from = 0,to = max(set_size(comb_mat)),by = set.size.by)),
                annotation_name_offset = unit(1,"cm")),
                row_names_gp = gpar(fontsize = 10),
                )

    draw_upset <- draw(upset_plot)
    order <- column_order(draw_upset)
    size <- comb_size(comb_mat,)

    decorate_annotation("intersection_size", {
        grid.text(size[order],
                  x = seq(size),
                  y = unit(size[order],
                           "native") +
                      unit(2, "pt"),
                  gp = gpar(fontsize = 10,
                            fontface = "plain"),
                  just = "bottom",
                  default.units = "native")
    })

    grid.text(label = label,gp = gpar(fontsize = 18),y = unit(label.vertical.position,"cm"),
              just = "top")

    dev.off()
}
