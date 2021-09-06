#' @title ModuleIC
#' @description Weighted gene co-expression network analysis for describing the correlation patterns among genes
#' @param vec.of.modulecolors a vector of module colors
#' @param my.moduleColors a module color object from WGCNA
#' @param my.IC a interconnectivity object from WGCNA
#' @param my.ExpData a dataset, features (genes, proteins ect) as columns and samples as rows.
#' @param my.n a Number indicating top n% most interconnected features in dataset
#' @param my.name a name of output plot(s)
#' @param my.softPower ??
#' @export
#' @import WGCNA
#' @seealso
#' @return
#' @examples \dontrun{
#' ...
#' }
#'

ModuleIC <- function(vec.of.modulecolors, my.moduleColors, my.IC, my.ExpData, my.n, my.softPower, my.name) {
    my.modules <- list()
    for (idx in 1:length(vec.of.modulecolors)) {
        moduleC <- colnames(my.ExpData)[which(my.moduleColors == vec.of.modulecolors[[idx]])]
        if (length(moduleC) <= 100) {
            pdf(paste0(my.name, "_module", vec.of.modulecolors[[idx]], "_moduleHM.pdf"), height=14, width=14)
            plotNetworkHeatmap(my.ExpData, moduleC, weights = NULL, useTOM = TRUE, power = my.softPower, networkType = "unsigned", main = "Heatmap of the network")
            dev.off()
        }
        moduleCIC <- my.IC[rownames(my.IC) %in% moduleC, ]
        moduleCIC <- moduleCIC[order(moduleCIC$kWithin,decreasing = TRUE), ]
        modTOP <- ceiling((length(moduleC)/100) * my.n)
        modTOP <- moduleCIC[1:modTOP,]
        modTOP$Module <- paste0("module", rep(vec.of.modulecolors[[idx]],nrow(modTOP)))

        modTOP$gene <- as.factor(rownames(modTOP))
        bp <- ggplot(modTOP, aes(x=gene, y=kWithin))+ geom_bar(stat="identity", fill=as.character(vec.of.modulecolors[[idx]]), width = 0.4) + theme_minimal() + geom_text(aes(label=gene), color = "grey30", size = 2.5, angle = 90, hjust = -.05) + theme(axis.text.x=element_blank())
        ggsave(paste0(my.name,"_", vec.of.modulecolors[[idx]], "_moduleIC.pdf"), plot = bp, width = 14, height = 10)

        my.modules[[idx]] <- modTOP

    }
    names(my.modules) <- vec.of.modulecolors
    return(my.modules)
}
