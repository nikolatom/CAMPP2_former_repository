#' @title A function for LASSO/Elastic net/Ridge regression
#' @description The function runs LASSO/Elastic net/Ridge regression
#' (depending on a hyperparameter alpha value). Parameter Lambda is
#' automatically estimated using k-fold cross-validation (cv.glmnet) by default.
#' The function optionally calculates also double cross-validation using
#' a test set. aValidation (test) set represents 1/4 of the samples from each
#' group.
#' Results are provided in the form of:
#' 1) the ames of the signifficant features passing filters based on their
#' coefficient values
#' 2) classification error rate (calculated if "validation=TRUE")
#' 3) cv.glmnet object.
#' @param seed an integer vector containing the random number generator (RNG)
#' state. Default is 123.
#' @param data a data frame of expression/abundance counts.
#' It's recommended to use normalized and batch corrected data.
#' @param group a factor specifying group for each sample (e.g. could be
#' represented by a column from a metadata file)
#' @param alpha a numeric vector specifying hyperparameter alpha. This value
#' must be set to 0.0 < x < 1.0 for Elastic Net (0.5 is default) or to 1.0 for
#' LASSO regression or to 0.0 for Ridge regression.
#' @param validation a boolean value for performing double cross-validation.
#' Validation (test) set represents 1/4 of the samples in each group.
#' As a result, cross validation error rate is provided. Default value is FALSE.
#' @param min.coef a numeric vector specifying a threshold for features'
#' filtering (e.g. genes) based on the coefficients which are calculated during
#' model fitting. Default value is > 0.
#' @param nfolds a numeric vector describing number of folds during Lambda
#' estimation which is based on a cross-validation. Although nfolds can be
#' as large as the sample size (leave-one-out CV), it is not recommended for
#' large datasets. Smallest value allowable is nfolds=3. Default is 10.
#' @export
#' @import glmnet
#' @import parallel
#' @import doMC
#' @import rngtools
#' @return a list of:
#' 1) coef.ma.names - a matrix of feature names having coefficients of best model
#' passing the filters (threshold defined by min.coef)
#' 2) class.error - logical/numeric value describing miss classification error
#' rate
#' 3) fit - cv.glmnet object
#' 4) coef.ma - a matrix of the features' coefficients (best model)
#' passing the filters (threshold defined by min.coef)
#' @examples {
#' LASSOFeature(seed=123, data=campp2_brca_1,
#' group=campp2_brca_1_meta$diagnosis, alpha=0.5, validation=TRUE,
#' min.coef = 0, nfolds=10)
#' }


    LASSOFeature <- function(seed=123, data, group, alpha=0.5, validation=FALSE, min.coef=0, nfolds = 10) {

        #estimate "family" parameter for glmnet
        if(length(levels(as.factor(group)))>2){
            family="multinomial"
        }else{
            family="binomial"
        }
        print(paste0(family," model selected for cv.glmnet."))


        if (validation == TRUE) {  ###calculating double cross-validation

            ll <- list()
            llev <- levels(as.factor(group))

            ##creates a list of indexes for each group
            for (idx in 1:length(llev)) {
                pos <- which(group == as.character(llev[idx]))
                ll[[idx]] <- pos
            }

            ##create a random sample selection (1/4 of the samples from each group in the input data). This is used only for the calculation of classification error rates
            samp <- unlist(lapply(ll, function(x) sample(x, ceiling((length(x)/4)))))

            data.val <- data[,samp]
            group.val <- as.character(group[samp])

            data <- data[,-samp]
            group <- as.character(group[-samp])
        }

        if(family == "multinomial") {
            RNGseed(seed)
            nth <- detectCores(logical = TRUE)
            registerDoMC(cores=nth)
            fit <- cv.glmnet(x = t(data), y = group, family="multinomial", type.multinomial = "grouped", nfolds = nfolds, alpha = alpha, parallel=TRUE) # An alternative described here: https://glmnet.stanford.edu/articles/glmnet.html
            # fit <- glmnet(x = t(data), y = group, family="multinomial", type.multinomial = "grouped", lambda = fit.cv$lambda.min, nfolds = nfolds, alpha = alpha, parallel=TRUE) #An alternative which is commonly used: https://www.statology.org/lasso-regression-in-r/; https://www.pluralsight.com/guides/linear-lasso-and-ridge-regression-with-r; https://www.r-bloggers.com/2020/05/quick-tutorial-on-lasso-regression-with-example/amp/
            coef <- coef(fit, s = "lambda.min")
            # coef <- coef(fit)
            coef.ma <- as(coef[[1]], "matrix")

        } else if (family == "binomial") {
            RNGseed(seed)
            fit <- cv.glmnet(x = t(data), y = group, family = "binomial", type.measure = "class", nfolds = nfolds, alpha = alpha)
            # fit <- glmnet(x = t(data), y = group, alpha = alpha, lambda = fit.cv$lambda.min, family = "binomial", type.measure = "class")

            coef <- coef(fit, s = "lambda.min")
            # coef <- coef(fit)
            coef.ma <- as(coef, "matrix")
        }

        ###Validation on test data
        if (validation == TRUE) {
            ## Here we calculate miss classification error; the fit from the data is applied on data.val
            class.error <- round(as.numeric(mean(predict(fit, t(data.val), s="lambda.min", type="class") != group.val))*100, digits = 2)
            cat(paste0("\nClassification error rate calculated on test data is = ", class.error, "%\n"))

        } else {
            class.error=NA
            cat("\nClassification error rate on test data was not selected.\n")
        }


        coef.ma.names <- names(coef.ma[coef.ma[,1] >min.coef, ]) ##original - obtain only names
        coef.ma <- as.data.frame(coef.ma[coef.ma[,1] >min.coef, ])[,1] #coefficients

        ## return genes and coefficients, classification error rates and model
        return(list("coef.ma.names"=coef.ma.names, "class.error"=class.error, "fit"=fit, "coef.ma"=coef.ma))

    }


