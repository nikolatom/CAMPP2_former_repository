#' @title Normalization, Filtering and Transformation
#' @description Function for Normalization, Filtering and Transformation. This pipeline does not handle background correction of single-channel intensity data or within array normalization two-color intensity data. See limma manual section on array normalization for more on this. Data may be fully normalized with limma (R/Rstudio) or another software and the pipeline re-run.
#' @param data a dataframe with expression/abundance counts
#' @param group a dataframe with samples groups
#' @param data.original a original data
#' @param sdata a dataframe with expression/abundance counts for second dataset
#' @param sgroup a dataframe with samples groups for second dataset
#' @param sdata.original a original data for second dataset
#' @param transform a transformation type
#' @param standardize a standardization type
#' @param technology a technology
#' @export
#' @seealso
#' @return normalized, filtered and transformed data
#' @examples \dontrun{
#' ...
#' }


applyNormalization<-function(data,sdata=NULL,data.original,sdata.original=NULL,group,sgroup=NULL,transform,standardize,technology){

    if (exists("data.original")) {
        data <- NormalizeData(technology[1], data, group, transform[1], standardize[1], data.original)
    } else {
        data <- NormalizeData(technology[1], data, group, transform[1], standardize[1])
    }


    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Second Dataset

    if(!is.null(sdata)) {
        if (length(technology) < 2) {
            stop("\nTwo datasets are input for correlation analysis, BUT argument technology only has length one. Length of technology must be two.\n")
        }
        if (length(transform) < 2) {
            stop("\nTwo datasets are input for correlation analysis, BUT argument transform only has length one. Length of transformmust be two, see.\n")
        }
    }



    if (!is.null(sdata)) {
        if (exists("sdata.original")) {
            sdata <- NormalizeData(technology[2], sdata, sgroup, transform[2], standardize[2], sdata.original)
        } else {
            sdata <- NormalizeData(technology[2], sdata, sgroup, transform[2], standardize[2])
        }
    }
    my_list<-list("data"=data,"sdata"=sdata)
    return(my_list)
}
