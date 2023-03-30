#' @title A wrapper function for random forest classification and variable selection
#' @description The function runs random forest classification and variable
#' selection by running the RFApply function. The classification and variable
#' selection are calculated 10 times for 10 random seeds and the intersection of
#' the selected variables (e.g., genes) is obtained.
#' The function optionally performs validation using a validation
#' dataset of the random forest classification if there is enough samples in the
#' dataset. This is determined by the split.size argument. If validation is not
#' performed, only variable selection using random forest is carried out.
#' The function also intersects the selected features across all 10 seed runs
#' as found during the feature selection process.
#' @param data a matrix of (transformed and normalized) feature counts from
#' "seq", "array", "ms" or "other" technology (with feature IDs as row names
#' and sample IDs as columns). If available, batch corrected data should be used.
#' @param group a factor specifying group for each sample (e.g. could be
#' represented by a column from a metadata file)
#' @param split.size an integer specifying the minimum number of samples that
#' the groups must contain in order to carry out random forest classification
#' and validation
#' @param test.train.ratio a floating point number between 0 and 1 representing
#' the ratio of samples to keep as validation dataset. For example, a
#' test.train.ratio = 0.25 splits 25 percent of the data into a validation dataset,
#' meaning 75 percent of the data will be kept as the training dataset.
#' @param num.trees.init an integer specifying number of trees to use for the first
#' forest in the feature selection process. Default is 5000.
#' @param num.trees.iterat an integer specifying number of trees to use for
#' all additional forests in the feature selection process. Default is 2000.
#' @export
#' @import caret
#' @return a list containing two elements:
#' 1) The first element is a list containing the output from the RFApply function
#' 2) The second element is a character vector containing the intersection of
#' selected variables (e.g., genes) from the feature selection process
#' @examples \dontrun{
#' campp2_brca_1_batchCorrected<-BatchCorrect(data=campp2_brca_1_normalized,
#' batch=campp2_brca_1_meta$tumor_stage,group=campp2_brca_1_meta$diagnosis,
#' technology="seq")
#'
#' campp2_brca_1_run_rf <-
#' RunRF(data = campp2_brca_1_batchCorrected,
#' group = campp2_brca_1_meta$diagnosis,
#' split.size = 5,
#' test.train.ratio = 0.25,
#' num.trees.init = 5000,
#' num.trees.iterat = 2000)
#' }

RunRF <- function(data,
                  group,
                  split.size,
                  test.train.ratio,
                  num.trees.init = 5000,
                  num.trees.iterat = 2000) {

    ## Run random forest

    # Length of each group
    group.size <- as.numeric(table(group))
    test.train <- unique(group.size >= split.size)

    if (FALSE %in% test.train) {

        # Run random forest variable selection
        print("Validation is not going to be done due to the low number (<split.size) of samples in at least one of the sample groups. Only variable selection will be performed.")
        RF.results <- RFApply(data = data, group = group, validation = FALSE, test.train.ratio = test.train.ratio,
                              num.trees.init = num.trees.init, num.trees.iterat = num.trees.iterat)

    } else {

        # Run random forest classification and variable selection
        RF.results <- RFApply(data = data, group = group, validation = TRUE, test.train.ratio = test.train.ratio,
                              num.trees.init = num.trees.init, num.trees.iterat = num.trees.iterat)

    }

    # Intersect features selected in each seed run
    intersect.sel.vars <- Reduce(intersect, RF.results$Sel.vars)

    return(list("RFResults" = RF.results, "VarsSelect" = intersect.sel.vars))

}



