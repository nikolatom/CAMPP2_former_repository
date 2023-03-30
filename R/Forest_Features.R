#' @title Random forest fitting and variable selection
#' @description This function implements random forest on feature counts
#' (e.g. genes) and allows for feature selection using the random forest
#' algorithm. For the feature selection process, the recommended value for
#' number of trees is 5000 for the first forest and 2000 for all additional
#' forests which are used as default values in this function.
#' In the variable selection process, variables are eliminated iteratively
#' by excluding the least important variables from each random forest. 20 percent of
#' the variables are excluded following each iteration. The out-of-bag (OOB)
#' error is used as criterion for determining the final selected variables.
#' Besides variable selection, a random forest model is also fitted which
#' is used for classification of the samples.
#' The input data is split into training and validation datasets where the
#' training data is used for variable selection and classification and the
#' random forest classifier is validated on the validation dataset.
#' Splitting of the data is determined by the test.train.ratio. For example,
#' a test.train.ratio = 0.25 splits 25 percent of the data into a validation dataset,
#' meaning 75 percent of the data will be kept as the training dataset.
#' @param seed an integer containing a random seed number. Default is 123.
#' @param data a matrix of (transformed and normalized) feature counts from
#' "seq", "array", "ms" or "other" technology (with feature IDs as row names
#' and sample IDs as columns). If available, batch corrected data should be used.
#' @param group a factor specifying group for each sample (e.g. could be
#' represented by a column from a metadata file)
#' @param validation a boolean indicating if validation will be performed
#' on test data. If TRUE validation will be performed on test data. If FALSE
#' validation will not be performed on test data. Default is FALSE.
#' @param test.train.ratio a floating point number between 0 and 1 representing
#' the ratio of samples to keep as validation dataset. For example, a
#' test.train.ratio = 0.25 splits 25 percent of the data into a validation dataset,
#' meaning 75 percent of the data will be kept as the training dataset.
#' @param num.trees.init an integer specifying number of trees to use for the first
#' forest in the feature selection process. Default is 5000.
#' @param num.trees.iterat an integer specifying number of trees to use for
#' all additional forests in the feature selection process. Default is 2000.
#' @export
#' @import varSelRF
#' @importFrom randomForest randomForest predict
#' @import caret
#' @import rngtools
#' @return a list of four elements: 1) a varSelRF object containing results of variable
#' selection using random forest, 2) a randomForest object containing random forest model
#' fitted to data, 3) a factor containing predictions of test data using fitted random
#' forest model and 4) a confusionMatrix object containing confusion matrix of test data
#' using fitted random forest model.
#' If validation = FALSE, the last three elements in output will be NA.
#' @examples \dontrun{
#' campp2_brca_1_batchCorrected<-BatchCorrect(data=campp2_brca_1_normalized,
#' batch=campp2_brca_1_meta$tumor_stage,group=campp2_brca_1_meta$diagnosis,
#' technology="seq")
#'
#' campp2_brca_1_forest_features <-
#' ForestFeatures(seed = 123,
#' data = campp2_brca_1_batchCorrected,
#' group = campp2_brca_1_meta$diagnosis,
#' validation = TRUE,
#' test.train.ratio = 0.25,
#' num.trees.init = 5000,
#' num.trees.iterat = 2000)
#' }

ForestFeatures <- function(seed = 123,
                           data,
                           group,
                           validation = FALSE,
                           test.train.ratio,
                           num.trees.init = 5000,
                           num.trees.iterat = 2000) {

    group <- as.factor(group)

    if (validation == TRUE) {

        ll <- list()
        llev <- levels(as.factor(group))

        # Generate a list of length number of groups
        # The first element in the list is the indexes of the samples belonging to group 1
        # The n-th element in the list is the indexes of the samples belonging to group n
        # where n is the number of groups
        for (idx in 1:length(llev)) {

            pos <- which(group == as.character(llev[idx]))
            ll[[idx]] <- pos

        }

        # Split randomly a subset of the samples into a validation and training set
        # A ratio of test.train.ratio is used for the splitting
        # These samples will be used as validation data set
        # It will contain samples from both groups
        RNGseed(seed)
        samp <- unlist(lapply(ll,
                              function(x) sample(x,
                                                 ceiling((length(x) * test.train.ratio)))))
        # Generate the validation data set
        data.val <- t(data[,samp])
        group.val <- group[samp]

        # Generate the training data set, i.e. the original data minus the validation data
        data.train <- t(data[,-samp])
        group.train <- group[-samp]

        # Carry out variable selection using random forest on training data and out-of-bag errors
        RNGseed(seed)
        RFvars <- varSelRF(data.train,
                           group.train,
                           ntree = num.trees.init,
                           ntreeIterat = num.trees.iterat,
                           vars.drop.frac = 0.2)

        # Fit random forest algorithm on training data
        RNGseed(seed)
        RFforest <- randomForest(data.train,
                                 group.train,
                                 ntree = num.trees.init,
                                 ntreeIterat = num.trees.iterat,
                                 vars.drop.frac = 0.2)

        # Use random forest to predict test data
        pred <- predict(RFforest,
                        newdata = data.val)

        # Create confusion matrix using predicted group classes and actual group classes on test data
        RFpred <- confusionMatrix(data = pred,
                                  reference = group.val)

        # Combine results of variable selection, predictions, confusion matrix and random forest
        # into one list
        res <- list(RF.var.sel = RFvars,
                    RF.model = RFforest,
                    RF.pred = pred,
                    RF.confusion.mat = RFpred)

    } else {

        # Transpose data
        data <- t(data)

        # Carry out variable selection using random forest and out-of-bag errors
        RNGseed(seed)
        RFvars <- varSelRF(data,
                           group,
                           ntree = num.trees.init,
                           ntreeIterat = num.trees.iterat,
                           vars.drop.frac = 0.2)

        # Return results from variable selection using random forest
        res <- list(RFvars,
                    NA,
                    NA,
                    NA)
    }

    return(res)

}







