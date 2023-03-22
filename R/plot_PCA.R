#' @title Principal Component Analysis (PCA)
#' @description A function calculating PCA. PCA plot shows the relationships
#' between samples (based on squared euclidean distances) in the data set.
#' Each dot represents one sample. This function provides also a scree plot
#' describing percentages of explained variances by each PCA component and
#' plots with contributions of variables to PCA1 and PCA2.
#' @param data a dataframe of expression/abundance counts
#' @param group a factor specifying group for each sample (e.g. could be
#' represented by a column from a metadata file)
#' @param cols a vector of colors (one color for each group)
#' @param show.PCA.labels a boolean value (TRUE or FALSE) specifying if elements
#' (e.g. samples) should be labelled in the PCA plot. Labeling is based on
#' column names of the input data. Default value is FALSE.
#' @param prefix a character defining a prefix of output file
#' @param scale a boolean, if TRUE then data are scaled to unit variance. Default is FALSE
#' @export
#' @import factoextra
#' @import FactoMineR
#' @return 1) scree plot; 2) plot of contributions of variables to PC1; 3) plot
#' of contributions of variables to PC2; 4) 2D PCA plot projecting samples over
#' first 2 principal components
#' @examples {
#' PCAPlot(campp2_brca_1_batchCorrected, as.factor(campp2_brca_1_meta$subtype), show.PCA.labels =FALSE, cols=NULL, prefix="test_PCA")
#' }

PCAPlot <- function(data, group, show.PCA.labels =FALSE, cols=NULL, prefix="", scale=FALSE) {


        ###parse TRUE/FALSE into "all"/"none".
        if(show.PCA.labels==TRUE){
            show.PCA.labels<-"all"
        } else if (show.PCA.labels==FALSE) {
            show.PCA.labels=="none"
        } else {
            stop(paste0("The value ", show.PCA.labels, " defined as show.PCA.labels parameter is not supported. Supported values are TRUE/FALSE."))
        }

        res.pca <- PCA(t(data),  graph = FALSE, ncp=10, scale = scale) # principal component analysis
        fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50), ggtheme = theme_classic(), ncp=20) # visualize eigenvalues/variances
        ggsave(paste0(prefix,"_screeplot.png"))

        fviz_contrib(res.pca, choice = "var", axes = 1, top = 15, ggtheme = theme_classic()) # contributions of variables to PC1
        ggsave(paste0(prefix,"_contributions_PC1.png"))

        fviz_contrib(res.pca, choice = "var", axes = 2, top = 15, ggtheme = theme_classic()) # contributions of variables to PC2
        ggsave(paste0(prefix,"_contributions_PC2.png"))

        # create PCA plot, visualize samples by groups
        fviz_pca_ind(res.pca,
                     label = show.PCA.labels, # hide individual labels
                     habillage = as.factor(group), # color by groups
                     palette = cols,
                     addEllipses = TRUE, # Concentration ellipses
                     repel=TRUE,
                     ggtheme = theme_classic(),
                     labelsize = 2
        )
        ggsave(paste0(prefix,"_PCA.png"))

}
