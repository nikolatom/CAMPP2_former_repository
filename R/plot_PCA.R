#' @title Principal Component Analysis (PCA)
#' @description A function calculating PCA. PCA plot shows the relationships
#' between samples (based on squared euclidean distances) in the data set.
#' Each dot represents one sample. This function provides also a scree plot
#' describing percentages of explained variances by each PCA component and
#' plots with contributions of variables to PCA1 and PCA2. Top 3 PCA components
#' are visualized in 3D plot.
#' @param data a dataframe of expression/abundance counts
#' @param group a factor specifying group for each sample (e.g. could be
#' represented by a column from a metadata file)
#' @param cols a vector of colors (one color for each group)
#' @param PCA.labels a text ("all", "none) specifying the elements to be
#' labelled. Default value is "none"
#' @param prefix a character defining a prefix of output file
#' @export
#' @import factoextra
#' @import FactoMineR
#' @import pca3d
#' @seealso
#' @return 1) scree plot, 2) contributions of variables to PC1, 3) contributions
#' of variables to PC2, 4) 2D PCA plot, 5) 3D PCA plot
#' @examples \dontrun{
#' PCAPlot(campp2_brca_1_batchCorrected, as.factor(campp2_brca_1_meta$subtype), PCA.labels ="none", cols=NULL, prefix="test_PCA")
#' }

PCAPlot <- function(data, group, PCA.labels ="none", cols=NULL, prefix="") {


        res.pca <- PCA(t(data),  graph = FALSE, ncp=10, scale = FALSE) # principal component analysis
        fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50), ggtheme = theme_classic(), ncp=20) # visualize eigenvalues/variances
        ggsave(paste0(prefix,"_screeplot.png"))

        fviz_contrib(res.pca, choice = "var", axes = 1, top = 15, ggtheme = theme_classic()) # contributions of variables to PC1
        ggsave(paste0(prefix,"_contributions_PC1.png"))

        fviz_contrib(res.pca, choice = "var", axes = 2, top = 15, ggtheme = theme_classic()) # contributions of variables to PC2
        ggsave(paste0(prefix,"_contributions_PC2.png"))

        # create PCA plot, visualize samples by groups
        fviz_pca_ind(res.pca,
                     label = PCA.labels, # hide individual labels
                     habillage = as.factor(group), # color by groups
                     palette = cols,
                     addEllipses = TRUE, # Concentration ellipses
                     repel=TRUE,
                     ggtheme = theme_classic(),
                     labelsize = 2
        )
        ggsave(paste0(prefix,"_PCA.png"))

        # PCA 3D
        pca<-prcomp(t(data), scale.=FALSE)
        pca3d(pca, group=as.factor(group), fancy=FALSE,legend="topleft",show.ellipses=TRUE, ellipse.ci=0.75, show.plane=FALSE)
        snapshotPCA3d(file=paste0(prefix,"_3D_PCA.png"))
}
