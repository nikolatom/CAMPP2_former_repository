#' @title Correlation analysis
#' @description If two datasets are provided (for correlation and/or network analysis), this option should be specified as a comma separated list (without quotes or parenthesis!) of length two,
#' first entry referring to data file 1 and second entry referring to the data file 2.
#' @param my.d1 first dataframe with values to be correlated. These must have the same dimensions and rownames must the same.
#' @param my.d2 second dataframe with values to be correlated. These must have the same dimensions and rownames must the same.
#' @param my.filename a name of output plot
#' @export
#' @import ggplot2
#' @import plyr
#' @seealso
#' @return correlation plots
#' @examples \dontrun{
#' ...
#' }

number_ticks <- function(n) {function(limits) pretty(limits, n)}

CorrAnalysis <- function(d1, d2, my.filename) {

    my.names <- rownames(d1)

    d1 <- as.matrix(d1)
    d2 <- as.matrix(d2)

    pear_corr <- sapply(1:nrow(d1), function(i) cor(d1[i,], d2[i,], method = "pearson"))

    # pearson correlation p-values
    pear_p_val <- sapply(1:nrow(d1), function(i) cor.test(d1[i,], d2[i,], method = "pearson")$p.value)

    # correction for multiple testing with fdr
    fdr <- p.adjust(pear_p_val, method = "fdr")

    # make dataframe
    pear_corr_full <- data.frame(my.names, pear_corr, pear_p_val, fdr, log2(1/fdr))
    colnames(pear_corr_full) <- c("name", "cor_coef", "pval", "fdr", "Inverse_Scaled_FDR")

    corr.plot <- ggplot(pear_corr_full, aes(name, cor_coef)) +  geom_point(aes(colour = Inverse_Scaled_FDR), size=7, stroke = 0, shape = 16) + scale_color_viridis(begin = 0.2, end = 0.8, direction = -1, option="cividis") + scale_y_continuous(breaks=number_ticks(10)) + theme_bw() +  theme(panel.grid.major.x = element_blank()) + geom_text(aes(label=name), size=4, hjust = 0.8, vjust=-0.2, color="grey30") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y = element_text(size=14, color = "black"), axis.title = element_text(size=16, color = "black"), legend.text = element_text(size=16), legend.title = element_text(size=14)) + xlab("") + ylab("Correlation Coefficient") + geom_hline(yintercept=c(0.0, 0.5, -0.5), color=rep("grey30", 3))
    ggsave(paste0(my.filename, "_corrplot.pdf"), bg = "transparent", width=12, height=6, plot = corr.plot)
    return(pear_corr_full)
}
