#' @title Estimate K-means clusters
#' @description A function to estimate the number of K-means.
#' Number of K-means is based on the Bayesian information criterion (BIC)
#' provided by mclust package. A summary describing the best model is printed on
#' the screen during the calculation.
#' @param data a data frame of feature (e.g. gene) counts
#' @export
#' @import mclust
#' @return a list including:
#' 1) a data frame with a number of clusters corresponding to G from BIC object.
#' G represents a number of mixture components in the model corresponding to
#' the optimal BIC
#' 2) an 'mclustBIC' object, which is the result of applying mclustBIC to data
#' @examples {
#' EstimateKmeans(t(campp2_brca_1_batchCorrected[1:3000,]))
#' }

EstimateKmeans <- function(data) {
    BIC <- mclustBIC(data)
    mod1 <- Mclust(data, x = BIC) #obtained model
    num.km.clusters <- as.numeric(mod1$G)  #number of clusters according to the best model
    print(paste0("number of clusters according to the best model, G: ",num.km.clusters))
    summary.clust<-summary(mod1, parameters = FALSE) #summary for the best model
    print(summary.clust)
    return(list("num.km.clusters"=num.km.clusters, "BIC"=BIC))
}
