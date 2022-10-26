#' @title Plotting Distributions
#' @description Plotting the counts' Distributions
#' @param data a data.frame of expression/abundance counts. N.B only a subset
#' of variables should be input, not intended for the full expression matrix!
#' (See function "FitDistributions")
#' @param fitted.data a list of lists obtained as output from the function
#' "FitDistributions". The data and fitted.data must have the same dimensions,
#' e.g. length(fitted.data) == nrows(data).
#' @export
#' @import ggplot2
#' @import squash
#' @import viridis
#' @seealso
#' @return distribution plots
#' @examples \dontrun{
#'
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

