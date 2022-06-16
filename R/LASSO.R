#' @title LASSO function
#' @description Run LASSO function to obtain DA features
#' @param my.seed a randomly generated vector of seeds
#' @param my.data a data marix of values
#' @param my.group a vector of integers specifying group
#' @param my.LAorEN a hyperparameter alpha. This value must be set to 0.0 < x < 1.0 for Elastic Net (0.5 is default) or to 1.0 for LASSO regression.
#' @param my.validation perform validation TRUE/FALSE
#' @param my.multinorm perfom multinorm TRUE/FALSE; IF my.multinorm=TRUE, then analysis will be multinomial, IF my.multinorm=FALSE, then then analysis will be binomial
#' @export
#' @import glmnet
#' @seealso
#' @return DA/DE features from LASSO
#' @examples \dontrun{
#' ...
#' }


LASSOFeature <- function(my.seed, my.data, my.group, my.LAorEN, my.validation=FALSE, my.multinorm=TRUE) {

    if (my.validation == TRUE) {

        ll <- list()
        llev <- levels(as.factor(my.group))

        for (idx in 1:length(llev)) {
            pos <- which(my.group == as.character(llev[idx]))
            ll[[idx]] <- pos
        }

        my.samp <- unlist(lapply(ll, function(x) sample(x, ceiling((length(x)/4)))))


        my.data.val <- my.data[,my.samp]
        my.group.val <- as.integer(my.group[my.samp])

        my.data <- my.data[,-my.samp]
        my.group <- as.integer(my.group[-my.samp])
    }

    if(my.multinorm == TRUE) {
        set.seed(my.seed)
        nth <- detectCores(logical = TRUE)
        registerDoMC(cores=nth)
        my.fit <- cv.glmnet(x = t(my.data), y = my.group, family="multinomial", type.multinomial = "grouped", nfolds = 10, alpha = my.LAorEN, parallel=TRUE)
        my.coef <- coef(my.fit, s=my.fit$lambda.1se)
        my.ma <- as(my.coef[[1]], "matrix")
    } else {
        set.seed(my.seed)
        my.fit <- cv.glmnet(x = t(my.data), y = my.group, family = "binomial", type.measure = "class", nfolds = 10, alpha = my.LAorEN)
        my.coef <- coef(my.fit, s=my.fit$lambda.1se)
        my.ma <- as(my.coef, "matrix")
    }


    if (my.validation == TRUE) {
        meanerror <- round(as.numeric(mean(predict(my.fit, t(my.data.val), s=my.fit$lambda.1se, type="class") != my.group.val))*100, digits = 2)
        cat(paste0("\nOne LASSO/EN run out completed. Mean class error was = ", meanerror, "%\n"))

    } else {
        meanerror <- round(as.numeric(mean(predict(my.fit, t(my.data), s=my.fit$lambda.1se, type="class") != my.group))*100, digits = 2)
        cat(paste0("\nOne LASSO/EN run out completed. Mean cross-validation error was = ", meanerror, "%\n"))
    }


    my.ma <- names(my.ma[my.ma[,1] != 0, ])

    return(list(my.ma, meanerror))

    rm(my.fit)
    rm(my.coef)
}


