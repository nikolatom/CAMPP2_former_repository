#' @title Apply Normalization and Transformation
#' @description A function for applying Normalization ("mean" or "median") and Transformation ("log2", "logit", "log10" and voom). This pipeline does not handle background correction of single-channel intensity data or within array normalization two-color intensity data. See limma manual section on array normalization for more on this. Data may be fully normalized with limma or another software and the pipeline re-run.
#' @param data a dataframe of expression/abundance counts
#' @param group a factor specifying sample group from metadata
#' @param data.original a original data (before removal of zeros)
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
    
    # Check if data contains zeros and negative values.
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
