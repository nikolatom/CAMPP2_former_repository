#' @title Normalization and Transformation
#' @description Function for Normalization ("mean" or "median") and Transformation ("log2", "logit" or "log10"). This pipeline does not handle background correction of single-channel intensity data1 or within array normalization two-color intensity data. See limma manual section on array normalization for more on this. Data may be fully normalized with limma or another software and the pipeline re-run.
#' @param data1 a dataframe with expression/abundance counts (genes as rows; samples as columns)
#' @param group1 a dataframe with groups metadata for each sample (given from metadata table)
#' @param group2 a dataframe with groups metadata for each sample from a second dataset
#' @param data1.original a original data1 including zero values
#' @param data2 a dataframe with expression/abundance counts for a second dataset
#' @param data2.original a original data2 including zero values
#' @param transform a string vector of length 1 (or two in case of 2 datasets) defining transformation type for each dataset ("log2", "logit" or "log10").
#' @param standardize a string vector of length 1 (or two in case of 2 datasets) defining standardization method  ("mean" or "median"); By default, data from "seq" technology are normalized by "TMM"; "array" technology is standardized using "quantile"
#' @param technology a string vector of length 1 (or two in case of 2 datasets) defining technology used for generating the data. Allowed types are: 'array', 'seq', 'ms' or 'other'.
#' @export
#' @seealso
#' @return Elist of normalized, filtered and transformed (gene) counts data
#' @examples \dontrun{
#' ...
#' }sample(nrow(data$E), 10)



applyNormalization2<-function(data1,data2=NULL,data1.original=NULL,data2.original=NULL,group1,group2=NULL,transform,standardize,technology){

    data1 <- NormalizeData(data=data1, data.original=data1.original, group=group1, transform[1], standardize[1], technology[1])

    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Second Dataset

    if(!is.null(data2)) {
        if (length(technology) != 2 || length(technology) == 1) {
            stop("\nTwo datasets are defined in the analysis, BUT argument technology only has length one. Technology must be defined as string vector of length two.\n")
        }
        if (length(transform) != 2 || length(transform) == 1) {
            stop("\nTwo datasets are defined in the analysis, BUT argument transform only has length one. Transform must be defined as string vector of length two.\n")
        }


        data2 <- NormalizeData(data=data2, data.original=data2.original, group=group2, transform[2], standardize[2], technology[2])

    }
    return(list("data1"=data1,"data2"=data2))
}
