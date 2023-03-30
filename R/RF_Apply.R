#' @title Applying random forest across seeds
#' @description This function applies the ForestFeatures function using 10 random
#' seeds. The ForestFeatures function performs random forest on feature counts
#' (e.g. genes) and allows for feature selection using the random forest
#' algorithm. This function extracts various results from the feature selection
#' process such as the selected variables and the OOB errors. Moreover, the OOB
#' error which resulted in the final selected model and the final selected
#' variables is extracted. The final selected model is the one with the smallest
#' number of genes whose error rate is within 1 standard error of the minimum
#' error rate of all forests. The function also extracts results from the process
#' of random forest fitting to the input data such as OOB error, predictions on
#' test data and corresponding accuracies. Finally, the function calculates the
#' average initial importance of each variable across the seed runs.
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
#' @import randomForest
#' @import caret
#' @return a list of seven elements: 1) a list where each element is a character string
#' containing selected variables from the random forest feature selection process from
#' each seed run, 2) a list where each element is out-of-bag errors (of type numeric)
#' from all random forest models used in the feature selection process for each seed run,
#' 3) a list where each element is the OOB error (of type numeric) from the best random
#' forest model (i.e. the final selected model) of the feature selection process for each
#' seed run, 4) a matrix containing average initial importance (before any variable deletion)
#' for each feature across seed runs, 5) a list where each element is the OOB error (of type
#' numeric) from fitted random forest model for each seed run (if a random forest was not
#' fitted then this will be NA), 6) a matrix ccontaining accuracy and 95 percent confidence
#' interval for predictions of test data using fitted random forest model where each row
#' in the matrix is a seed run, and 7) a vector of integers containing seeds used in the
#' procedure.
#' @examples \dontrun{
#' campp2_brca_1_batchCorrected<-BatchCorrect(data=campp2_brca_1_normalized,
#' batch=campp2_brca_1_meta$tumor_stage,group=campp2_brca_1_meta$diagnosis,
#' technology="seq")
#'
#' campp2_brca_1_rf_apply <-
#' RFApply(data = campp2_brca_1_batchCorrected,
#' group = campp2_brca_1_meta$diagnosis,
#' validation = TRUE,
#' test.train.ratio = 0.25,
#' num.trees.init = 5000,
#' num.trees.iterat = 2000)
#' }

RFApply <- function(data,
                    group,
                    validation = FALSE,
                    test.train.ratio,
                    num.trees.init = 5000,
                    num.trees.iterat = 2000) {

    # Get randomly 10 seeds
    seeds <- sample(1:1000, 10)

    # Assign empty lists to store results
    RFvars <- list() # store selected variables from feature selection process
    RFvars.oob <- list() # store out-of-bag errors from feature selection process
    RFsel.vars.oob <- list() # store out-of-bag error from the best random forest model (i.e. the final selected model)
    RFimp <- list() # store initial importance of variables before any variable deletion from feature selection process
    RFoob.fit <- list() # store OOB error from fitted random forest model
    RFacc <- list() # store accuracy and 95% confidence interval for predictions of test data using fitted random forest model

    # For each seed
    for (idx in 1:length(seeds)) {

        # Run random forest
        RF <- ForestFeatures(seed = seeds[[idx]],
                             data = data,
                             group = group,
                             validation = validation,
                             test.train.ratio = test.train.ratio,
                             num.trees.init = num.trees.init,
                             num.trees.iterat = num.trees.iterat)

        ## Extract different results from the random forest variable selection process

        # Obtain selected variables found from running feature selection using random forest
        RFvars[[idx]] <- RF[[1]]$selected.vars

        # Obtain OOB error from all random forest models used in feature selection process
        RFvars.oob[[idx]] <- RF[[1]]$selec.history$OOB

        # Obtain OOB error from the best random forest model (i.e. the final selected model)
        tmp <- which(RF[[1]]$selec.history$Number.Variables ==
                         RF[[1]]$best.model.nvars) # Find index of the RF that is the final selected model
        RFsel.vars.oob[[idx]] <- RF[[1]]$selec.history$OOB[tmp] # Get OOB error of final selected RF model

        # Obtain initial importance of features before any variable deletion from feature selection process
        RFimp[[idx]] <- RF[[1]]$initialImportances


        ## Extract results from the process of fitting a random forest model to the data

        # If a random forest model was fitted
        if (is(class(RF[[2]]), "randomForest")){
            # Obtain OOB error from fitted random forest model
            RFoob.fit[[idx]] <- tail(RF[[2]]$err.rate, n = 1)[,1]

        } else {

            RFoob.fit[[idx]] <- NA

        }

        # If predictions were performed on test data
        if (is(class(RF[[4]]), "confusionMatrix")){
            # Obtain accuracy and 95% confidence interval for predictions of test data using fitted random forest model
            RFacc[[idx]] <- RF[[4]]$overall[c(1, 3, 4)]

        } else {

            RFacc[[idx]] <- NA

        }

        # Print message stating how many seed runs are completed
        print(paste0(idx, " out of ", length(seeds), " runs completed."))

    }

    ## Calculate average importance for each feature across seed runs

    # Initialize all.df to be importance of all features from first seed run
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
                  RFvars.oob,
                  RFsel.vars.oob,
                  RFimp,
                  RFoob.fit,
                  do.call(rbind,
                          RFacc),
                  seeds)

    # Assign names to elements in list
    names(RFres) <- c("Sel.vars",
                      "Var.sel.oob",
                      "Sel.rf.oob",
                      "Mean.importance.features",
                      "oob.rf.model",
                      "accuracy.rf.model",
                      "seeds")

    return(RFres)

}









