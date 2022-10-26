#' @title Random forest fitting and variable selection
#' @description This function implements random forest on feature counts
#' (e.g. genes) and allows for feature selection using the random forest
#' algorithm.
#' @param seed an integer containing a random seed number.
#' @param data a matrix of (transformed and normalized) feature counts from
#' "seq", "array", "ms" or "other" technology (with feature IDs as row names
#' and sample IDs as columns).
#' @param group a vector of type factor containing the group ID for each
#' sample
#' @param validation a boolean indicating if validation will be performed
#' on test data. If TRUE validation will be performed on test data. If FALSE
#' validation will not be performed on test data. Defaults to FALSE.
#' @param num_trees an integer specifying number of trees to use for the first
#' forest in the feature selection process
#' @param num_trees_iterat an integer specifying number of trees to use for
#' all additional forests in the feature selection process
#' @export
#' @import varSelRF
#' @import randomForest
#' @import caret
#' @seealso
#' @return a list of four elements: 1) results of variable selection using random forest,
#' 2) random forest model fitted to data, 3) predictions of test data using fitted random
#' forest model and 4) confusion matrix of test data using fitted random forest model.
#' If validation = FALSE, the last three elements in output will be NA.
#' @examples \dontrun{
#' campp2_brca_1_forest_features <-
#' Forest_Features(seed = 173,
#' data = campp2_brca_1,
#' group = campp2_brca_1_meta$diagnosis,
#' validation = TRUE,
#' num_trees = 30,
#' num_trees_iterat = 15)
#' }

Forest_Features <- function(seed,
                            data,
                            group,
                            validation = FALSE,
                            num_trees,
                            num_trees_iterat) {

  group <- as.factor(group)

  if (validation == TRUE) {

    ll <- list()
    llev <- levels(as.factor(group))

    # Generate a list of length number of groups
    # The first element in the list is the indexes of the samples belonging to group 1
    # The second element in the list is the indexes of the samples belonging to group 2
    for (idx in 1:length(llev)) {

      pos <- which(group == as.character(llev[idx]))
      ll[[idx]] <- pos

    }

    # Generate randomly a vector of 1/4 of the number of samples
    # These samples will be used as validation data set
    # It will contain samples from both groups
    set.seed(5)
    samp <- unlist(lapply(ll,
                             function(x) sample(x,
                                                ceiling((length(x) / 4)))))
    # Generate the validation data set
    data.val <- t(data[,samp])
    group.val <- group[samp]

    # Generate the training data set, i.e. the original data minus the validation data
    data <- t(data[,-samp])
    group <- group[-samp]

    # Carry out variable selection using random forest on training data and out-of-bag errors
    set.seed(seed)
    RFvars <- varSelRF(data,
                       group,
                       ntree = num_trees,
                       ntreeIterat = num_trees_iterat,
                       vars.drop.frac = 0.2)

    # Fit random forest algorithm on training data
    set.seed(seed)
    RFforest <- randomForest(data,
                             group,
                             ntree = num_trees,
                             ntreeIterat = num_trees_iterat,
                             vars.drop.frac = 0.2)

    # Use random forest to predict test data
    pred <- predict(RFforest,
                    newdata = data.val)

    # Create confusion matrix using predicted group classes and actual group classes on test data
    RFpred <- confusionMatrix(data = pred,
                              reference = group.val)

    # Combine results of variable selection, predictions, confusion matrix and random forest
    # into one list
    res <- list(RF_var_sel = RFvars,
                RF_model = RFforest,
                RF_pred = pred,
                RF_confusion_mat = RFpred)

  } else {

    # Transpose data
    data <- t(data)

    # Carry out variable selection using random forest and out-of-bag errors
    set.seed(seed)
    RFvars <- varSelRF(data,
                       group,
                       ntree = num_trees,
                       ntreeIterat = num_trees_iterat,
                       vars.drop.frac = 0.2)

    # Return results from variable selection using random forest
    res <- list(RFvars,
                NA,
                NA,
                NA)
  }

  return(res)

}







