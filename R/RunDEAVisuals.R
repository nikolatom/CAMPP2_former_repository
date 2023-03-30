#' @title A wrapper for a visualisation of differential expression/abundance
#' analysis
#' @description This function wraps several functions for the visualisation
#' of DEA results:
#' 1) A function for making volcano plots from DEA results.
#' The input is typically an output from RunDEA (re-formatted matrix of
#' differential expression/abundance results from limma, please, check RunDEA
#' for more details). By default, top 15 up- and down-regulated features
#' will be labeled on the Volcano plot. Important: Features are selected from
#' the table with all the groups' comparisons.
#' 2) A function for making Venn diagrams to visualize the size of
#' intersections of the features between the groups
#' (e.g., differentially expressed genes for each subtype)
#' 3) A function for making Upset plots to visualize the size of
#' intersections of the features present each group
#' (e.g., differentially expressed genes for each subtype)
#' @param data A data frame from RunDEA containing genes and their
#' corresponding AveExpr, logFC, p.value, adj.p.value ect. Mandatory columns
#' are "logFC", "adj.P.Val" and "HUGO_ID." HUGO_ID could be obtained using
#' AddGeneName function.
#' @param prefix a prefix for the output filename.
#' @param cutoff.FDR a numeric value for adj.P.Val cutoff. Default = 0.01
#' @param cutoff.logFC a numeric value for logFC cutoff. Default = 1
#' @param n.labeled.features a number of the most up-/down- regulated features
#' labeled in the volcano plot. Default = 15
#' @param control.group A string vector defining control group name for each
#' dataset (e.g., "healthy" or "normal").
#' @export
#' @import ggplot2
#' @import ggrepel
#' @import VennDiagram
#' @import viridis
#' @import ComplexHeatmap
#' @return several .png files with:
#' 1) Volcano plot of differentially expressed features saved into .png
#' 2) Venn Diagram showcasing the size of intersections of the features
#' between the sample groups
#' 3) Upset plot showcasing the size of intersections of the features
#' between the sample groups
#' @examples {
#' RunDEAVisuals(campp2_brca_1_DEA_HUGO, cutoff.FDR = 0.05, cutoff.logFC = 1,
#' n.labeled.features = 15, control.group= "healthy", prefix="test_DEA_visuals")
#' }

RunDEAVisuals <- function(data, cutoff.FDR = 0.01, cutoff.logFC = 1, n.labeled.features = 15, control.group, prefix){

    print(" - Generating Volcano Plot.")
    MakeVolcano(data, prefix, cutoff.logFC, cutoff.FDR)

    ##Groups' analysis
    if (length(unique(data$comparison)) > 1){
        groups.full.list <- split.data.frame(data, data$comparison) ##data is already groups.full.list vs control group (based on grepl before)

        groups.feature.list=list()

        for (features in groups.full.list){
            groups.feature.list <- append(groups.feature.list,list(features$name))  ## Get feature names (ENSEMBL IDs) for each comparison
        }

        names(groups.feature.list) <- gsub("-","",gsub(control.group,"",names(groups.full.list))) ##Get groups.full.list without control group and use it as names of feature lists


        print(" - Printing Venn Diagram.")
        MakeVennDiagram(prefix,groups.feature.list)

        print(" - Printing Upset Plot.")
        MakeUpset(prefix,groups.feature.list)

        ##Subtypes and direction-aware expression
        groups.full.list.updown<-split.data.frame(data,list(data$comparison,data$dir))

        groups.feature.list.updown=list()

        for (features in groups.full.list.updown){
            groups.feature.list.updown <- append(groups.feature.list.updown,list(features$name))
        }

        names(groups.feature.list.updown) <- gsub("-","",gsub(control.group,"",names(groups.full.list.updown))) ##Get groups.full.list without control group and use it as names of feature lists

        print(" - Printing direction-wise Upset Plot.")
        MakeUpset(paste0(prefix,".updown"),groups.feature.list.updown)
        # print(" - Printing Venn Diagram for the first dataset.")
        # MakeVennDiagram(paste0(prefix,".updown"),groups.feature.list.updown)  ##this is not working due to one argument not used by format 'Incorrect number of elements.'
    }
}
