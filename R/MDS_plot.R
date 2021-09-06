#' @title Multidimensional Scaling Plot
#' @description Data overview. The plot shows the relationship between samples (squared euclidian distance) in the dataset. Each dot represents one sample.
#' @param my.data a dataframe of expression/abundance counts
#' @param my.group a vector of integers specifying the group (may be the same or different, length should be equal to ncol(dataframe))
#' @param my.lables a vector of IDs for labeling (may be the same or different, length should be equal to ncol(dataframe))
#' @param my.cols a vector of colors (one color for each group)
#' @export
#' @import heatmap.plus
#' @import ggplot2
#' @import squash
#' @import viridis
#' @seealso
#' @return distribution plots
#' @examples \dontrun{
#' MDSPlot (counts, groups, labels, cols)
#' }
#
# library("RColorBrewer")
# my.cols<- brewer.pal(4,"Set1")
# my.group<-as.factor(metadata[,10])
# my.labels<-metadata[,1]
# my.data<-counts2
# mds_plot<-MDSPlot(counts2,my.group, my.labels, my.cols)
# mds_plot

MDSPlot <- function(my.data, my.group, my.labels, my.cols) {
    d<-dist(t(my.data))
    fit <- cmdscale(d,eig=TRUE, k=2)
    res <-data.frame(names=rownames(fit$points),M1=fit$points[,1],M2=fit$points[,2])
    ggplot(res, aes(x=M1, y=M2)) + geom_point(aes(fill = my.group, colour = my.group), shape=21, size = 5, stroke = 0.1) + geom_text(aes(x=M1,y=M2, label= my.labels)) + guides(fill = guide_legend(override.aes = list(shape = 22))) + scale_fill_manual(values=my.cols) + scale_colour_manual(values=rep("white", length(my.cols))) + theme_bw() + theme(legend.title=element_blank()) + theme(legend.text = element_text(size = 16, face="bold"), axis.title=element_text(size=16,face="bold")) + theme(legend.position = "top") + theme(axis.text=element_text(size=16, face="bold")) + theme(axis.text = element_text(colour = "black"))
}

