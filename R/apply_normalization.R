#' @title Apply Normalization and Transformation
#' @description A function for applying Normalization ("mean" or "median") and Transformation ("log2", "logit", "log10" and voom). This pipeline does not handle background correction of single-channel intensity data1 or within array normalization two-color intensity data. See limma manual section on array normalization for more on this. Data may be fully normalized with limma or another software and the pipeline re-run.
#' @param data1 a dataframe with expression/abundance counts (genes as rows; samples as columns)
#' @param group1 a dataframe with groups metadata for each sample (given from metadata table)
#' @param group2 a dataframe with groups metadata for each sample from a second dataset
#' @param data1.original a original data1 including zero values
#' @param data2 a dataframe with expression/abundance counts for a second dataset
#' @param data2.original a original data2 including zero values
#' @param transform a string vector of length 1 (or two in case of 2 datasets) defining transformation type for each dataset ("log2", "logit", "log10"). voom transform is automatically applied on "seq" data.
#' @param standardize a string vector of length 1 (or two in case of 2 datasets) defining standardization method  ("mean" or "median"); By default, data from "seq" technology are normalized by "TMM"; "array" technology is standardized using "quantile"
#' @param technology a string vector of length 1 (or two in case of 2 datasets) defining technology used for generating the data. Allowed types are: 'array', 'seq', 'ms' or 'other'.
#' @export
#' @seealso
#' @return Elist of normalized, filtered and transformed (gene) counts data
#' @examples \dontrun{
#' ...
#' }
#'


applyNormalization<-function(data,data.original=NULL,group,transform,standardize,technology){

    # Check if data1 contains zeros and negative values.
    hasZero <- unique(as.vector(data == 0))
    hasNeg <- unique(as.vector(data < 0))
    if (TRUE %in% hasNeg) {
        print("data includes negative values and cannot be log transformed.")
    }
    if (TRUE %in% hasZero) {
        print("data include zero values which might be replaced by running function ReplaceZerosApply")
    }

 
    ###Run normalization

    data <- NormalizeData(data, data.original, group, transform, standardize, technology)

    return(data)
}
