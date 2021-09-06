#' @title Estimate Kmeans
#' @description Function for Estimating Kmeans; Number of kmeans will be based on the bayesian information criterion(BIC); RNAseq with
#' many genes, multiple samples of 3000 variables will be generated
#' and tested to overcome issues with computational time and the
#' consensus of best n kmeans will be returned.
#' Clustering may only be performed one dataset at a time!
#' @param my.data a dataframe of expression/abundance counts
#' @param n = number of sample subsets to generate (The number of clusters tested will be based on number of samples, fewer samples will result in fewer kmeans tested.)
#' @export
#' @import mclust
#' @seealso
#' @return clustering results and plots
#' @examples \dontrun{
#' EstimateKmeans(counts, n)
#' }
#
# library("mclust")
# df<-counts2
# n<-3
# kmeans<-EstimateKmeans(df,n)

EstimateKmeans <- function(df, n) {
    BIC <- mclustBIC(df)
    mod1 <- Mclust(df, x = BIC)
    Ks <- as.numeric(summary(mod1, parameters = TRUE)$G)
    cat(paste0("\nCluster run complete - out of ", length(n), " in total..."))
    return(Ks)
}
