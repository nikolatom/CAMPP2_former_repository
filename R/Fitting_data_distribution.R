#' @title Fitting Data Distributions
#' @description A function for fitting data distribution
#' @param data a data frame of expression/abundance counts, N.B only a subset of
#' variables should be input, not intended for the full expression matrix!
#' Depending on the character of the data, various distributions are checked.
#' In case of only integer numbers present in the dataset, poison and normal
#' distributions are tested. In case of presence of float values in the dataset,
#' weibull, gamma, log normal and normal distributions are tested. In case of
#' negative float values, only normal distribution is calculated.
#' @export
#' @import fitdistrplus
#' @seealso
#' @return a list of the results from fitdist function
#' @examples \dontrun{
#' campp2_brca_1_distributions <- FitDistributions(campp2_brca_1[1:10,])
#' }


FitDistributions <- function(data) {
    discretetype <- unique(as.numeric(apply(data, 1, function(x) x%%1==0)))  ## check for integers
    hasNeg <- unique(as.numeric(data < 0))  ##check for negative values
    list.of.lists <- list()
    for(idx in 1:nrow(data)) {
        if (FALSE %in% discretetype) {
            if (TRUE %in% hasNeg) {
                try(fit_n <- fitdist(as.numeric(data[idx,]), "norm"), silent = TRUE)
                l <- list()
                if(exists("fit_n")) {
                    l[[length(l)+1]] <- fit_n
                }
                list.of.lists[[idx]] <-  l
            } else {
                try(fit_w  <- fitdist(as.numeric(data[idx,]), "weibull"), silent = TRUE)
                try(fit_g  <- fitdist(as.numeric(data[idx,]), "gamma"), silent = TRUE)
                try(fit_ln <- fitdist(as.numeric(data[idx,]), "lnorm"), silent = TRUE)
                try(fit_n <- fitdist(as.numeric(data[idx,]), "norm"), silent = TRUE)
                l <- list()
                if(exists("fit_w")) {
                    l[[length(l)+1]] <- fit_w
                }
                if(exists("fit_g")) {
                    l[[length(l)+1]] <- fit_g
                }
                if(exists("fit_ln")) {
                    l[[length(l)+1]] <- fit_ln
                }
                if(exists("fit_n")) {
                    l[[length(l)+1]] <- fit_n
                }
                list.of.lists[[idx]] <-  l
            }
        }
        if (discretetype == TRUE) {
            try(fit_p <- fitdist(as.numeric(data[idx,]), "pois"), silent = TRUE)
            try(fit_n <- fitdist(as.numeric(data[idx,]), "norm"), silent = TRUE)
            l <- list()
            if(exists("fit_p")) {
                l[[length(l)+1]] <- fit_p
            }
            if(exists("fit_n")) {
                l[[length(l)+1]] <- fit_n
            }
            list.of.lists[[idx]] <-  l
        }
    }
    names(list.of.lists) <- rownames(data)
    return(list.of.lists)
}
