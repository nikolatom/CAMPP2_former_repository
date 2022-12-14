#' @title A wrapper function for random forest classification and variable selection
#' @description The function runs random forest classification and variable
#' selection by running the RFApply function. The classification and variable
#' selection are calculated 10 times for 10 random seeds and the intersection of
#' the selected variables (e.g., genes) is obtained.
#' The function optionally performs validation using a validation
#' dataset of the random forest classification if there is enough samples in the
#' dataset. If validation is not performed, only variable selection using random
#' forest is carried out.
#' The function also intersects the selected features across all 10 seed runs
#' as found during the feature selection process.
#' @param data a matrix of (transformed and normalized) feature counts from
#' "seq", "array", "ms" or "other" technology (with feature IDs as row names
#' and sample IDs as columns)
#' @param group a factor specifying group for each sample (e.g. could be
#' represented by a column from a metadata file)
#' @param split.size an integer specifying the minimum number of samples that
#' the groups must contain in order to carry out random forest classification
#' and validation
#' @param num.trees.init an integer specifying number of trees to use for the first
#' forest in the feature selection process
#' @param num.trees.iterat an integer specifying number of trees to use for
#' all additional forests in the feature selection process
#' @export
#' @import varSelRF
#' @import randomForest
#' @import caret
#' @seealso
#' @return a list containing two elements:
#' 1) The first element is a list containing the output from the RFApply function
#' 2) The second element is a character vector containing the intersection of
#' selected variables (e.g., genes) from the feature selection process
#' @examples \dontrun{
#' campp2_brca_1_run_rf <-
#' RunRF(data = campp2_brca_1_batchCorrected,
#' group = campp2_brca_1_meta$diagnosis,
#' split.size = 5,
#' num.trees.init = 30,
#' num.trees.iterat = 15)
#' }

RunRF <- function(data, group, split.size, num.trees.init, num.trees.iterat) {

    ## Run random forest

    # Length of each group
    group.size <- as.numeric(table(group))
    test.train <- unique(group.size >= split.size)

    if (FALSE %in% test.train) {

        # Run random forest variable selection
        print("Validation is not going to be done due to the low number (<split.size) of samples in at least one of the sample groups. Only variable selection will be performed.")
        RF.results <- RFApply(data = data, group = group, validation = FALSE,
                              num.trees.init = num.trees.init, num.trees.iterat = num.trees.iterat)

    } else {

        # Run random forest classification and variable selection
        RF.results <- RFApply(data = data, group = group, validation = TRUE,
                              num.trees.init = num.trees.init, num.trees.iterat = num.trees.iterat)

    }

    # Intersect features selected in each seed run
    intersect.sel.vars <- Reduce(intersect, RF.results$Sel.vars)

    return(list("RFResults" = RF.results, "VarsSelect" = intersect.sel.vars))

}



