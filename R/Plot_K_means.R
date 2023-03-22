#' @title Plot Kmeans
#' @description Function for Plotting Kmeans clustering; A folder with MDS plots will be returned for the best n cluster, based on BIC
#' @param my.data a dataframe of expression/abundance counts
#' @param clus.list a list of k-means?
#' @param k A number of clusters tested will be based on number of samples, fewer samples will result in fewer kmeans tested.
#' @param my.labels metadata column for sample labeling
#' @param my.name ?dataset name?
#' @export
#' @import mclust
#' @import rngtools
#' @seealso
#' @return clustering results and plots
#' @examples \dontrun{
#' ...
#' }

PlotKmeans <- function(my.data, clus.list, k, my.labels, my.name) {
    res.list <- list()
    nclus <- unique(unlist(clus.list)) ##doesn't make sense as it generates only "1"
#    nclus <- length(unlist(clus.list))

    if (unique(is.na(nclus)) == TRUE) {
        cat("\nNo 'best' ks could be determined. There may be little or poor clustering of samples. Ks 1:5 will be returned.\n")
        nclus <- k ###MIGHT BE GOOD TO SETUP THIS CONDITION NCLUS!=1 (or smth similar)
    }
    nclus <- sort(nclus)
    for (idx in 1:length(nclus)) {
        RNGseed(10)
#        nth <- detectCores(logical = TRUE)
        Kclus <- kmeans(t(my.data), nclus[[idx]])
        Clusters <- as.factor(paste0("C",data.frame(Kclus$cluster)$Kclus.cluster))
        d<-dist(t(my.data))
        fit <- cmdscale(d,eig=TRUE, k=2)
        res <-data.frame(names=rownames(fit$points),M1=fit$points[,1],M2=fit$points[,2])
        res$Clusters <- Clusters
        res.list[[idx]] <- Clusters
        Pal <- viridisLite::viridis(length(levels(as.factor(res$Clusters))))
        p <- ggplot(res, aes(x=M1, y=M2)) + geom_point(aes(fill = Clusters, colour = Clusters), shape=21, size = 3, stroke = 0.1) + guides(fill = guide_legend(override.aes = list(shape = 22))) + scale_fill_manual(values=Pal) + scale_colour_manual(values=rep("white", length(Pal))) + theme_bw() + geom_text(label = my.labels, colour = "grey50", nudge_x = 0.25, nudge_y = 0.25, size =3) + theme(legend.title=element_blank()) + theme(legend.text = element_text(size = 14), axis.title=element_text(size=14)) + theme(legend.position = "top") + theme(axis.text=element_text(size=14)) + theme(axis.text = element_text(colour = "black"))
        ggsave(file=paste0("BestKmeans_C", as.character(nclus[[idx]]), ".pdf"), p, width = 10, height = 8)
    }
    names(res.list) <- paste0("Clus", nclus)
    res.df <- as.data.frame(res.list)
    return(res.df)
}
