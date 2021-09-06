#' @title WGCNA Analysis
#' @description Weighted gene co-expression network analysis for describing the correlation patterns among genes
#' @param my.data a dataframe with expression/abundance values
#' @param my.thresholds a vector of length two. First argument specifies the minimum size of a module and the second argument specifies the dissimilarity cutoff for merging of modules.
#' @param my.name a name of output plot(s)
#' @export
#' @import WGCNA
#' @seealso
#' @return
#' @examples \dontrun{
#' ...
#' }

WGCNAAnalysis <- function(my.data, my.thresholds, my.name) {

    enableWGCNAThreads()
    powers <-  c(c(1:10), seq(from = 12, to=20, by=2));
    sft <- pickSoftThreshold(my.data, dataIsExpr = TRUE,powerVector = powers,corFnc = cor,corOptions = list(use = 'p'),networkType = "unsigned")

    pdf(paste0(my.name,"_WGCNA_softpowerplot.pdf"))
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n")
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],col="red")
    dev.off()

    # Set power to best softpower estimate
    softPower <- sft$powerEstimate

    if (is.na(softPower) | softPower > 9) {
        if (nrow(my.data) < 20) {
            softPower <- 9
        } else if (nrow(my.data) >= 20 & nrow(my.data) < 30) {
            softPower <- 8
        } else if (nrow(my.data) >= 30 & nrow(my.data) < 40) {
            softPower <- 7
        } else {
            softPower <- 6
        }
    }

    # Construct adjecancy and topology matrix
    # ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    # Adjacency matrix and Topology Matrix
    adj <-adjacency(my.data, power = softPower)

    if(ncol(my.data) > 8000) {

        cat("Dataset too large for classic WGCNA, running blockwiseModules intead.\n")

        WGCNAres <- blockwiseModules(
            my.data,
            weights = NULL,
            checkMissingData = TRUE,
            blocks = NULL,
            maxBlockSize = 8000,
            blockSizePenaltyPower = 5,
            nPreclusteringCenters = as.integer(min(ncol(my.data)/20, 100*ncol(my.data)/5000)),
            randomSeed = 12345,
            corType = "pearson",
            maxPOutliers = 1,
            quickCor = 0,
            pearsonFallback = "individual",
            power = softPower,
            networkType = "unsigned",
            replaceMissingAdjacencies = FALSE,
            suppressTOMForZeroAdjacencies = FALSE,
            TOMType = "unsigned",
            TOMDenom = "min",
            getTOMs = NULL,
            saveTOMs = FALSE,
            saveTOMFileBase = "blockwiseTOM",
            deepSplit = 2,
            detectCutHeight = 0.995,
            minModuleSize = my.thresholds[1],
            minBranchEigennodeDissim = mergeCutHeight,
            tabilityCriterion = c("Individual fraction", "Common fraction"),
            reassignThreshold = 1e-6,
            minCoreKME = 0.5,
            minCoreKMESize = my.thresholds[1]/3,
            minKMEtoStay = 0.3,
            mergeCutHeight = 0.25)

        moduleColors <- WGCNAres$colors

    } else {
        TOM <- TOMsimilarity(adj)
        colnames(TOM) <- colnames(my.data)
        rownames(TOM) <- colnames(my.data)

        cat("Done with topology calculations.\n")

        dissTOM <- (1-TOM)

        # Dendogram
        geneTree <- hclust(as.dist(dissTOM),method="ward.D2")

        cat("Done with hclust.\n")

        # Construct Modules
        # ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = my.thresholds[1])


        # Extract module color labels
        dynamicColors <- labels2colors(dynamicMods)

        # Change colors to viridis
        my.cols <- data.frame(dynamicColors)
        colnames(my.cols) <- c("oldcols")

        n.colors <- length(levels(as.factor(dynamicColors)))

        WGCNA.cols <- c("#5C6D70", "#E88873", "#B5BD89", "#E8AE68", "#104F55", "#4062BB", "#72195A", "#8ACB88", "#A997DF", "#4B296B", "#BC69AA", "#1C3144", "#C3D898", "#C2C1A5", "#B23A48", "#3A7D44", "#B6F9C9", "#FF8360", "#296EB4", "#E8C1C5", "#2E282A", "#AAA95A", "#34D1BF", "#3454D1", "#F2F4CB", "#898980", "#ADF5FF", "#372549", "#7DD181", "#F56476")

        if(n.colors > 30) {
            l <- n.colors-30
            WGCNA.cols <- c(WGCNA.cols, viridis(l, option="inferno"))
        }

        #colortransform <- data.frame(levels(as.factor(dynamicColors)), viridis(length(levels(as.factor(dynamicColors))), option="inferno"))
        colortransform <- data.frame(levels(as.factor(dynamicColors)), WGCNA.cols[1:n.colors])
        colnames(colortransform) <- c("oldcols", "colortransform")
        dynamicColors <- as.character(join(my.cols, colortransform)$colortransform)

        # Eigen features in each module
        MEList <- moduleEigengenes(my.data, colors = dynamicColors)
        MEs <- MEList$eigengenes
        MEDiss <- (1-cor(MEs))

        # Cluster module eigengenes
        if (ncol(MEDiss) <= 1) {
            print("N.B. You only have one module in your tree, please consider changing parameter -x. Specify a smaller minimum module size (default is 10) or decrease % dissimilarity for mergining modules (default is 25%).")
        } else {
            METree <- hclust(dist(MEDiss), method = "ward.D2")
        }

        cat("Done with hclust of MEDiss.\n")

        merge <- mergeCloseModules(my.data, dynamicColors, cutHeight = my.thresholds[2]/100, verbose = 3)
        mergedColors <- merge$colors
        mergedMEs <- merge$newMEs

        # Plot merged modules
        pdf(paste0(my.name,"_WGCNA_moduleTree.pdf"), height = 10, width = 12)
        plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),c("Dynamic Tree Cut", "Merged dynamic"),dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
        dev.off()

        moduleColors <- mergedColors
        colorOrder <- c("grey", standardColors(50))
        moduleLabels <- match(moduleColors, colorOrder)-1
    }


    # Interconnectivity of features in modules
    # ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    IC <- intramodularConnectivity(adj,  moduleColors)
    mod.cols <- levels(as.factor(moduleColors))

    cat("Done with IC.\n")

    WGCNAres <- ModuleIC(mod.cols, moduleColors, IC, my.data, my.thresholds[3], softPower, my.name)
    return(WGCNAres)
}
