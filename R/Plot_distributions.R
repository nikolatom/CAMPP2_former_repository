#' @title Plot Distributions
#' @description Plotting the feature counts' Distributions.
#' As results, multiple plots are obtained for each feature:
#' 1) Cullen and Frey graph
#' 2) Histogram and theoretical densities
#' 3) Empirical and theoretical CDFs
#' 4) Q-Q plot
#' 5) P-P plot
#'
#' For some of the figures, plotted distributions depends on the characteristics
#' of the input data:
#' In case of only integer numbers are present in the dataset, poison and normal
#' distributions are tested. In case of presence of float values,
#' Weibull, gamma, log normal and normal distributions are tested. In case of
#' negative float values, only normal distribution is calculated.
#'
#' Resulting figures are saved into the .pdf file.
#' @param data a data.frame of feature counts, the same which was
#' used as an input for "FitDistributions" function.
#' @param fitted.data a list of distribution descriptions obtained as an output
#' from the function "FitDistributions" where the same feature counts are used
#' as an input.
#' @export
#' @import ggplot2
#' @import squash
#' @import viridis
#' @return a .pdf file including multiple plots:
#' 1) Cullen and Frey graph
#' 2) Histogram and theoretical densities
#' 3) Empirical and theoretical CDFs
#' 4) Q-Q plot
#' 5) P-P plot
#' @examples {
#' PlotDistributions(campp2_brca_1_batchCorrected[1:10,], campp2_brca_1_distributionsFit)
#' }

PlotDistributions <- function(data, fitted.data) {
    discretetype <- unique(as.vector(apply(data, 1, function(x) x%%1==0)))
    hasNeg <- unique(as.vector(data < 0))
    for(idx in 1:length(fitted.data)) {
        pdf(paste0(names(fitted.data)[idx], ".pdf"), height = 8, width = 12)
        par(mfrow=c(2,3))
        if (FALSE %in% discretetype) {
            if (TRUE %in% hasNeg) {
                descdist(as.numeric(data[idx,]), discrete = FALSE, boot = 500, obs.col = viridis(1), boot.col = viridis(5)[4])
                plot.legend <- c("norm")
                plot.colors <- viridis(1)
            } else {
                descdist(as.numeric(data[idx,]), discrete = FALSE, boot = 500, obs.col = viridis(1), boot.col = viridis(5)[4])
                plot.legend <- c("Weibull", "lognormal", "gamma", "norm")

                plot.colors <- viridis(4)
            }
        }
        if (!FALSE %in% discretetype) {
            descdist(as.numeric(data[idx,]), discrete = TRUE,  boot = 500, obs.col = viridis(1), boot.col = viridis(5)[4])
            plot.legend <- c("poisson", "norm")
            plot.colors <- viridis(2)
        }
        denscomp(fitted.data[[idx]], legendtext = plot.legend, fitcol = plot.colors)
        cdfcomp (fitted.data[[idx]], legendtext = plot.legend, fitcol = plot.colors, datapch=16)
        qqcomp  (fitted.data[[idx]], legendtext = plot.legend, fitcol = plot.colors, fitpch=16)
        ppcomp  (fitted.data[[idx]], legendtext = plot.legend, fitcol = plot.colors, fitpch=16)
        dev.off()
    }
}

