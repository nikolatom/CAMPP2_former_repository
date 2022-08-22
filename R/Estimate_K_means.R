#' @title Estimate K-means clusters
#' @description A function to estimate the number of K-means.
#' Number of K-means is based on the Bayesian information criterion (BIC)
#' provided by mclust package. A summary describing the best model is printed on
#' the screen during the calculation.
#' @param data a data frame of feature (e.g. gene) counts
#' @export
#' @import mclust
#' @seealso
#' @return a data frame with a number of clusters
#' @examples \dontrun{
#' EstimateKmeans(t(campp2_brca_1_batchCorrected[1:2000,]))
#' }

EstimateKmeans <- function(data) {
    BIC <- mclustBIC(data)
    plot(BIC)
    mod1 <- Mclust(data, x = BIC) #obtained model
    # fviz_nbclust(my_data, kmeans, method = "silhouette") #alternative approach
    num.km.clusters <- as.numeric(mod1$G)  #number of clusters according to the best model
    print(paste0("number of clusters according to the best model, G: ",num.km.clusters,"\n"))
    summary.clust<-summary(mod1, parameters = FALSE) #summary for the best model
    print(summary.clust)
    return(num.km.clusters)

}
