#' @title K-means clustering
#' @description A wrapper function for K-means clustering.
#' For data sets with >1000 features, multiple subsets (number of sets = number
#' of features/1000; rounded up to the higher integer) will be generated,
#' maximum is 10 subsets. A number of randomly selected features (e.g. genes)
#' in 1 subset is limited to 2000.
#' A number of clusters for each subset is automatically estimated using BIC and
#' a consensus (based on all the subsets) of best n k-means is returned and used
#' for clustering.
#' If K-means estimation using BIC fails, a number of clusters will be based on
#' the number of samples (2-6 clusters for data set with less than 100 samples;
#' 2-11 for data sets with 101-500 samples, and 2-16 clusters for data sets with
#' more than 500 samples).
#' As a result, an information about the clusters assigned to each sample and
#' the details from PCA are provided.
#' The clusters are visualized on a PCA plot and saved into the png.
#' A summary describing the best model is printed on the screen during
#' calculation.
#' @param data a data.frame of feature (e.g. gene) counts
#' @param PCA.labels a string vector ("all" or "none") specifying which elements
#' (e.g. samples) should be labelled (labeling is based on column names of the
#' input data). Default value is "none".
#' @param cols a vector of colors (one color for each group)
#' @param prefix a character string defining a prefix of output file.
#' @param num.km.clusters a vector of manually defined number(s) of clusters.
#' By default, the values(s) are calculated automatically (a default value is
#' NULL).
#' @export
#' @import factoextra
#' @import FactoMineR
#' @seealso
#' @return a data.frame with cluster information assigned to each sample;
#' a list of results from PCA;
#' PCA plot(s) saved into png.
#' @examples \dontrun{
#' runKmeans(campp2_brca_1_batchCorrected, PCA.labels = "none", cols=NULL,
#' prefix="test", num.km.clusters=NULL)
#' }

runKmeans <- function(data, PCA.labels = FALSE, cols=NULL, prefix=NULL, num.km.clusters=NULL){

     if(is.null(prefix)){
        stop(print("Please, provide a prefix for the result files."))
     }

    if(is.null(num.km.clusters)){
        # Number of sample subsets to generate
        n.subsets <- 1:ceiling(nrow(data)/1000)
        if(length(n.subsets) > 10) {
            n.subsets <- 1:10
        }

        # Number of variables (genes) in each sample subset, limit is 2000.
        setsize <- nrow(data)
        if (setsize > 2000) {
            setsize <- 2000
        }

        # Number of kmeans to try if k-means estimation using BIC is not working
        if(ncol(data) <= 100) {
            nks <- 2:6
        } else if (ncol(data) > 100 && ncol(data) <= 500) {
            nks <- 2:11
        } else {
            nks <- 2:16
        }


        list.of.dfs <- list()

        for (idx in 1:length(n.subsets)) {
            df <- t(data[sample(nrow(data), setsize), ])  #make random gene selections
            list.of.dfs[[idx]] <- df  #create list of random selections
        }
        cat(paste0("\n", length(n.subsets), " clustering runs will be done in total...\n"))
        n.clusters.list <- lapply(list.of.dfs, function(x) EstimateKmeans(x)) #estimate how many clusters in the data subsets


        nclus <- unique(unlist(n.clusters.list))

        if (unique(is.na(nclus)) == TRUE) {
            cat(paste0("Number of clusters could not be determined using BIC. There may be little or poor clustering of samples. Alternatively, based on size of dataset, ", length(n.subsets), " sample sets will be generated of size ", setsize, " and ", length(nks), " clusters will be tested. \nRunning..."))

            nclus <- nks
        }
    }else{
        nclus <- num.km.clusters
    }

    nclus <- sort(nclus)
    paste0("Numbers of clusters being tested: ",nclus)

    res.pca <- PCA(t(data),  graph = FALSE, ncp=10, scale = FALSE) # principal component analysis
    res.list <- list()

    for (idx in 1:length(nclus)) {
        set.seed(10)
        Kclus <- kmeans(t(data), nclus[[idx]])
        Clusters <- as.factor(paste0("C",data.frame(Kclus$cluster)$Kclus.cluster))
        res.list[[idx]] <- Clusters
        if(length(cols) < nclus[[idx]]){
            print("Number of colours is defined by the number of groups and is smaller than the number of clusters. Colour scheme will be defined automatically.")
            cols <- NULL
        }
        fviz_pca_ind(res.pca,
                     label = PCA.labels, # show/hide individual labels; labels are taken from feature counts matrix automatically
                     habillage = as.factor(Clusters), # color by groups (clusters in this case)
                     palette = cols,
                     addEllipses = TRUE, # concentration ellipses
                     repel=TRUE,
                     ggtheme = theme_classic(),
                     title = paste0("k-means ", nclus[[idx]], " clusters"),
                     labelsize = 2

        )
        ggsave(paste0(prefix,"_PCA_Kmeans_C",nclus[[idx]],".png"))

    }
    names(res.list) <- paste0("Clus", nclus)
    res.clusters <- as.data.frame(res.list)

    return(list("res.clusters"=res.clusters,"res.pca"=res.pca))

}
