#' @title Multidimensional Scaling Plot
#' @description Data overview. The plot shows the relationship between samples (squared euclidian distance) in the dataset. Each dot represents one sample.
#' @param data a dataframe of expression/abundance counts
#' @param group a factor derived from metadata column selected as a sample group (e.g. diagnosis). Length should be equal to ncol(data)).
#' @param prefix A prefix for the output filename.
#' @param lables a vector of IDs for labeling (may be the same or different, length should be equal to ncol(dataframe))
#' @param cols a vector of colors (one color for each group)
#' @param MDS.labels a TRUE/FALSE statement specifying whether to include labels on the output plot. Default is TRUE (include labels).
#' @export
#' @import ggplot2
#' @import squash
#' @import viridis
#' @seealso
#' @return distribution plots
#' @examples \dontrun{
#' ...
#' }

MDSPlot <- function(data, group, prefix, labels, cols, MDS.labels=TRUE) {
    d<-dist(t(data))
    fit <- cmdscale(d,eig=TRUE, k=2)
    res <-data.frame(names=rownames(fit$points),M1=fit$points[,1],M2=fit$points[,2])
    if (!is.null(MDS.labels) & length(labels)==ncol(data)){
        mdsplot<-ggplot(res, aes(x=M1, y=M2)) + geom_point(aes(fill = group, colour = group),
                        shape=21, size = 5, stroke = 0.1) + geom_text(aes(x=M1,y=M2, label= labels)) +
                        guides(fill = guide_legend(override.aes = list(shape = 22))) + scale_fill_manual(values=cols) +
                        scale_colour_manual(values=rep("white", length(cols))) + theme_bw() + theme(legend.title=element_blank()) +
                        theme(legend.text = element_text(size = 16, face="bold"), axis.title=element_text(size=16,face="bold")) +
                        theme(legend.position = "top") + theme(axis.text=element_text(size=16, face="bold")) +
                        theme(axis.text = element_text(colour = "black"))
    } else if (length(labels)!=ncol(data)) {
        stop(" - The number of MDS labels does not match the number of samples in your input dataset. No MDS plot will be printed.")
    } else {
        mdsplot<-ggplot(res, aes(x=M1, y=M2)) + geom_point(aes(fill = group, colour = group), shape=21, size = 5, stroke = 0.1) +
                        guides(fill = guide_legend(override.aes = list(shape = 22))) + scale_fill_manual(values=cols) +
                        scale_colour_manual(values=rep("white", length(cols))) + theme_bw() + theme(legend.title=element_blank()) +
                        theme(legend.text = element_text(size = 16, face="bold"), axis.title=element_text(size=16,face="bold")) +
                        theme(legend.position = "top") + theme(axis.text=element_text(size=16, face="bold")) +
                        theme(axis.text = element_text(colour = "black"))
    }
    ggsave(paste0(prefix, "_MDSplot.png"), plot = mdsplot, dpi = 300, width = 8, height = 8)
}
