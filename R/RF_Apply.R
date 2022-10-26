#' @title Applying random forest across seeds
#' @description This function applies the Forest_Features function on 10 random
#' seeds. The Forest_Features function performs random forest on feature counts
#' (e.g. genes) and allows for feature selection using the random forest
#' algorithm.
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
#' @return a list of eight elements: 1) selected variables from random forest
#' feature selection process, 2) out-of-bag errors from random forest feature
#' selection process, 3) OOB error from first random forest model fitted during
#' feature selection process, 4) confusion matrix of predictions made based on
#' OOB data and first random forest model fitted during feature selection process,
#' 5) importance of features from first random forest model fitted during feature
#' selection process, 6) OOB error from fitted random forest model,
#' 7) accuracy and 95% confidence interval for predictions of test data using
#' fitted random forest model, and 8) seeds used in the procedure.
#' @examples \dontrun{
#' campp2_brca_1_rf_apply <-
#' RF_Apply(data = campp2_brca_1,
#' group = campp2_brca_1_meta$diagnosis,
#' validation = TRUE,
#' num_trees = 30,
#' num_trees_iterat = 15)
#' }

RF_Apply <- function(data,
                     group,
                     validation,
                     num_trees,
                     num_trees_iterat) {

  # Get randomly 10 seeds
  seeds <- sample(1:1000, 10)

  # Assign empty lists to store results
  RFvars <- list() # store selected variables from feature selection process
  RFvars_oob <- list() # store out-of-bag errors from feature selection process
  RFoob <- list() # store OOB error from first random forest model fitted during feature selection process
  RFclass <- list() # store confusion matrix of predictions made based on OOB data and first random forest model fitted during feature selection process
  RFimp <- list() # store importance of features from first random forest model fitted during feature selection process
  RFoob_fit <- list() # store OOB error from fitted random forest model
  RFacc <- list() # store accuracy and 95% confidence interval for predictions of test data using fitted random forest model

  # For each seed
  for (idx in 1:length(seeds)) {

    # Run random forest
    RF <- Forest_Features(seed = seeds[[idx]],
                          data = data,
                          group = group,
                          validation = validation,
                          num_trees = num_trees,
                          num_trees_iterat = num_trees_iterat)

    # Obtain selected variables found from running feature selection using random forest
    RFvars[[idx]] <- RF[[1]]$selected.vars

    # Obtain OOB error from all random forest models used in feature selection process
    RFvars_oob[[idx]] <- RF[[1]]$selec.history$OOB

    # Obtain median OOB error from first random forest model fitted during feature selection process
    RFoob[[idx]] <- median(RF[[1]]$firstForest$err.rate[,1])

    # Obtain confusion matrix of predictions made based on OOB data and first random forest model fitted during feature selection process
    RFclass[[idx]] <- RF[[1]]$firstForest$confusion

    # Obtain importance of features from first random forest model fitted during feature selection process
    RFimp[[idx]] <- RF[[1]]$firstForest$importance

    # If a random forest model was fitted
    if (class(RF[[2]]) == "randomForest") {

      # Obtain OOB error from fitted random forest model
      RFoob_fit[[idx]] <- tail(RF[[2]]$err.rate, n = 1)[,1]

    } else {

      RFoob_fit[[idx]] <- NA

    }

    # If predictions were performed on test data
    if (class(RF[[4]]) == "confusionMatrix") {

      # Obtain accuracy and 95% confidence interval for predictions of test data using fitted random forest model
      RFacc[[idx]] <- RF[[4]]$overall[c(1, 3, 4)]

    } else {

      RFacc[[idx]] <- NA

    }

    # Print message stating how many seed runs are completed
    print(paste0(idx, " out of ", length(seeds), " runs completed."))

  }

  ## Calculate average importance for each feature across seed runs

  # Initiatlize all.df to be importance of all features from first seed run
  all.df <- RFimp[[1]]

  # For each additional seed run
  for (idx in 2:length(RFimp)) {

    # Sum the importance of each feature across the seed runs
    all.df <- all.df + RFimp[[idx]]

  }

  # Divide the summed importance of each feature with the number of seed runs to get average
  RFimp <- all.df / length(RFimp)

  # Combine all results into one list
  RFres <- list(RFvars,
                RFvars_oob,
                RFoob,
                RFclass,
                RFimp,
                RFoob_fit,
                do.call(rbind,
                        RFacc),
                seeds)

  # Assign names to elements in list
  names(RFres) <- c("Sel_vars",
                    "Sel_vars_oob",
                    "Sel_vars_oob_first",
                    "Sel_vars_confusion_mat_first",
                    "Mean_importance_features_first",
                    "oob_rf_model",
                    "accuracy_rf_model",
                    "seeds")

  return(RFres)

}









