#' @title Correlation plots
#' @description Function for correlation scatter plots
#' @param my.data = a dataframe with expression/abundance counts for tissue or TIF
#' @param my.serumdata = a dataframe with expression/abundance counts for serum
#' @param my.filename = name of output plot
#' @param my.features ?
#' @export
#' @import ggplot2
#' @import plyr
#' @seealso
#' @return correlation plots
#' @examples \dontrun{
#' ...
#' }

CorrelationPlots <- function(my.data, my.serumdata, my.features, my.filename) {
    my.data <- my.data[rownames(my.data) %in% my.features,]
    my.serumdata <- my.serumdata[rownames(my.serumdata) %in% my.features,]
    features <- lapply(1:nrow(my.data), function(x) data.frame(as.numeric(my.data[x,]), as.numeric(my.serumdata[x,])))
    features <- lapply(features, setNames, c("TIF_Tissue", "Serum"))
    p1 <- lapply(features, function(x) ggplot(x, aes(TIF_Tissue, Serum)) + geom_point(shape=1, size=2.5) + theme_bw() + geom_smooth(method=lm, colour="black") +  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + theme(axis.text.x = element_text(size=12), axis.text.y = element_text(size=12, color = "black"), axis.title = element_text(size=16, color = "black"), legend.text = element_text(size=16)))
    names(p1) <- my.features
    for (i in 1:length(p1)) {
        p1[[i]] <- p1[[i]] + ggtitle(names(p1)[i])
    }
    nc <- ceiling(length(features)/2)
    pdf(paste0(my.filename,"_individual_corrplots.pdf"), height=6, width=12)
    multiplot(plotlist = p1, cols = nc)
    dev.off()
}
