#' @title Fitting Data Distributions
#' @description A function for fitting data distribution
#' @param my.data a dataframe of expression/abundance counts, N.B only a subset of variables should be input, not intended for the full expression matrix!
#' @export
#' @import fitdistrplus
#' @seealso
#' @return data for distribution plots
#' @examples \dontrun{
#' FitDistributions(counts)
#' }

# ##FOR TESTING PURPOSES ONLY
# library("fitdistrplus")
# counts4<-FitDistributions(counts2)
# #counts4<-FitDistributions(counts3)


FitDistributions <- function(my.data) {
    discretetype <- unique(as.vector(apply(my.data, 1, function(x) x%%1==0)))
    hasNeg <- unique(as.vector(my.data < 0))
    list.of.lists <- list()
    for(idx in 1:nrow(my.data)) {
        if (FALSE %in% discretetype) {
            if (TRUE %in% hasNeg) {
                try(fit_n <- fitdist(as.numeric(my.data[idx,]), "norm"), silent = TRUE)
                l <- list()
                if(exists("fit_n")) {
                    l[[length(l)+1]] <- fit_n
                }
                list.of.lists[[idx]] <-  l
            } else {
                try(fit_w  <- fitdist(as.numeric(my.data[idx,]), "weibull"), silent = TRUE)
                try(fit_g  <- fitdist(as.numeric(my.data[idx,]), "gamma"), silent = TRUE)
                try(fit_ln <- fitdist(as.numeric(my.data[idx,]), "lnorm"), silent = TRUE)
                try(fit_n <- fitdist(as.numeric(my.data[idx,]), "norm"), silent = TRUE)
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
            try(fit_p <- fitdist(as.numeric(my.data[idx,]), "pois"), silent = TRUE)
            try(fit_n <- fitdist(as.numeric(my.data[idx,]), "norm"), silent = TRUE)
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
    names(list.of.lists) <- rownames(my.data)
    return(list.of.lists)
}
