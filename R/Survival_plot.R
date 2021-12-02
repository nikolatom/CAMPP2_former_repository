#' @title Survival plot
#' @description Function for generating survival plot
#' @param my.survres a survival object in the form of a list with a cox regression for each feature.
#' @param my.filename a name of output plot
#' @param my.n ?
#' @export
#' @import survminer
#' @import survcomp
#' @seealso
#' @return survival plot
#' @examples \dontrun{
#' ...
#' }


SurvivalCOX <- function(my.survres, my.filename, my.n) {
    survival_conf  <- lapply(my.survres, function(x) summary(x)[2, c(4,6,7)])
    survival_conf  <- data.frame(do.call(rbind, survival_conf))
    colnames(survival_conf) <- c("HR", "lower", "upper")
    survival_pval <-  as.numeric(unlist(lapply(my.survres, function(x) anova(x)[1,3])))
    survival_conf$pval <- survival_pval
    survival_conf$fdr <- p.adjust(survival_pval, method = "fdr", n=length(survival_pval))
    survival_conf$feature <- as.factor(names(my.survres))
    survival_conf$Significant <- as.factor(ifelse(survival_conf$fdr <=0.05, 1, 0))
    survival_conf$InverseFDR <- 1/survival_conf$fdr


    if(nrow(survival_conf) > 100) {
        n.splits <- floor(nrow(survival_conf)/my.n)
        x <- 1:nrow(survival_conf)
        n <- n.splits
        n.splits <- split(x, cut(seq_along(x), n, labels = FALSE))
        features <- lapply(n.splits, function(x) survival_conf[x,])

        p1 <- lapply(features, function(x) ggplot(x, aes(feature, HR)) + geom_point(aes(colour = Significant, size = InverseFDR)) + scale_color_manual(values=c("grey70", "black")) + geom_errorbar(aes(ymax = upper, ymin = lower)) + geom_hline(yintercept=1) + theme_bw() + theme(axis.text.x = element_text(size=13, color = "black", angle = 90, hjust = 1), axis.text.y = element_text(size=12, color = "black"), axis.title = element_text(size=15)) + theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line( size=.1, color="black" )) + scale_y_continuous(breaks=pretty(limits,10)) + xlab("Features") + ylab("Hazard Ratio") + scale_y_continuous(trans = log2_trans(), breaks = trans_breaks("log2", function(x) 2^x),labels = trans_format("log2", math_format(2^.x))))
        for (i in 1:length(p1)) {
            pdf(paste0(as.character(names(p1)[i]),"_individual_corrplots.pdf"), height=6, width=12)
            print(p1[i])
            dev.off()
        }

    } else {
        pdf(paste0(my.filename, "_corrplot.pdf"), height=6, width=12)
        p1 <- ggplot(survival_conf, aes(feature, HR)) + geom_point(aes(colour = Significant, size = InverseFDR)) + scale_color_manual(values=c("grey70", "black")) + geom_errorbar(aes(ymax = upper, ymin = lower)) + geom_hline(yintercept=1) + theme_bw() + theme(axis.text.x = element_text(size=13, color = "black", angle = 90, hjust = 1), axis.text.y = element_text(size=12, color = "black"), axis.title = element_text(size=15)) + theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line( size=.1, color="black" )) + scale_y_continuous(breaks=number_ticks(10)) + xlab("Features") + ylab("Hazard Ratio") + scale_y_continuous(trans = log2_trans(), breaks = trans_breaks("log2", function(x) 2^x),labels = trans_format("log2", math_format(2^.x)))
        print(p1)
        dev.off()
    }
    return(survival_conf)
}
