#' @title Create model matrix
#' @description A function for creating a model matrix based on gene counts
#' metadata table and selected column names (covariates) from the metadata.
#' @param metadata Samples' metadata table should be imported using function
#' "import_metadata". Metadata must include exactly the same samples as gene
#' counts (data1) and samples must be sorted similarly.
#' @param group a factor specifying group for each sample (e.g. could be
#' represented by a column from a metadata file)
#' @param batch a factor specifying batch for each sample (e.g. could be
#' represented by a column from a metadata file). Default = NULL.
#' @param covariates a covariate(s) to include in the analysis provided
#' as a character vector. covarDEA can be defined only with a batch parameter
#' @export
#' @import limma
#' @import sva
#' @return a model matrix used e.g. for limma
#' @examples {
#' design.matrix<-Design.Matrix(campp2_brca_1_meta, campp2_brca_1_meta$diagnosis)
#' }

Design.Matrix <- function(metadata, group, batch=NULL, covariates=NULL) {

    if (is.null(batch) && is.null(covariates)) {
        design.matrix <- model.matrix(~0+group)
    } else if (!is.null(batch) && is.null(covariates)) {
        design.matrix <- model.matrix(~0+group+batch)
    } else if(!is.null(batch) && !is.null(covariates)) {
        names<-paste0(covariates, collapse = "+")
        design.matrix <- eval(parse(text=paste0("model.matrix(~0+group+batch","+",names,",data=metadata)")))
    }
    return(design.matrix)
}
