#' @title Plotting Distributions
#' @description Plotting the counts' Distributions
#' @param my.data a dataframe of expression/abundance counts, N.B only a subset of variables should be input, not intended for the full expression matrix!(See "Fitting Data Distributions" above)
#' @param list.of.lists an output from the function "Fitting Data Distributions". The my.data and list.of.list should have the same dimentions, e.g. length of list == nrows of dataframe.
#' @export
#' @import heatmap.plus
#' @import ggplot2
#' @import squash
#' @import viridis
#' @seealso
#' @return distribution plots
#' @examples \dontrun{
#' PlotDistributions (counts, listOfLists)
#' }


# ##FOR TESTING PURPOSES ONLY
# #during testing (without normalization) some of the genes/proteins were not ploted
#
# library("heatmap.plus")
# library("ggplot2")
# library("squash")
# library("heatmap.plus")
#
# #plot<-PlotDistributions(counts3, counts4) #not working with normalized data because Error in apply(my.data, 1, function(x) x%%1 == 0) :
# #dim(X) must have a positive length
#
# plot<-PlotDistributions(counts2, counts4)

PlotDistributions <- function(my.data, list.of.lists) {
    discretetype <- unique(as.vector(apply(my.data, 1, function(x) x%%1==0)))
    hasNeg <- unique(as.vector(my.data < 0))
    for(idx in 1:length(list.of.lists)) {
        pdf(paste0(names(list.of.lists)[idx], ".pdf"), height = 8, width = 12)
        par(mfrow=c(2,3))
        if (FALSE %in% discretetype) {
            if (TRUE %in% hasNeg) {
                descdist(as.numeric(my.data[idx,]), discrete = FALSE, boot = 500, obs.col = viridis(1), boot.col = viridis(5)[4])
                plot.legend <- c("norm")
                plot.colors <- viridis(1)
            } else {
                descdist(as.numeric(my.data[idx,]), discrete = FALSE, boot = 500, obs.col = viridis(1), boot.col = viridis(5)[4])
                plot.legend <- c("Weibull", "lognormal", "gamma", "norm")
                plot.colors <- viridis(4)
            }
        }
        if (!FALSE %in% discretetype) {
            descdist(as.vector(my.data[idx,]), discrete = TRUE,  boot = 500, obs.col = viridis(1), boot.col = viridis(5)[4])
            plot.legend <- c("poisson", "norm")
            plot.colors <- viridis(2)
        }
        denscomp(list.of.lists[[idx]], legendtext = plot.legend, fitcol = plot.colors)
        cdfcomp (list.of.lists[[idx]], legendtext = plot.legend, fitcol = plot.colors, datapch=16)
        qqcomp  (list.of.lists[[idx]], legendtext = plot.legend, fitcol = plot.colors, fitpch=16)
        ppcomp  (list.of.lists[[idx]], legendtext = plot.legend, fitcol = plot.colors, fitpch=16)
        dev.off()
    }
}

