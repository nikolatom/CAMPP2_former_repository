#' @title A wrapper for differential expression/abundance analysis
#' @description A function for running differential expression/abundance
#' analysis using a matrix of features (e.g., genes) counts as a mandatory
#' input. In the function, user is allowed to take advantage of defining multiple
#' covariates, thresholds for a log fold change (logFC), a false discovery
#' rate (FDR) and a blocking variable.
#' The differential expression/abundance analysis is based on voom
#' transformation and limma.
#' Results are provided in the form of tables with differentially expressed/
#' abundant features. Design a contrast matrices are also provided.
#' @param data a matrix of (transformed and normalized) feature counts from
#' "seq", "array", "ms" or "other" technology (with feature IDs as row names
#' and sample IDs as columns).
#' @param metadata a samples' metadata table should be imported using function
#' "import_metadata". Metadata must include exactly the same samples as gene
#' counts (data) and samples must be sorted similarly. In this function,
#' metadata are optional and are used to extract covariates for each sample.
#' Default = NULL.
#' @param batch a factor specifying batch for each sample (e.g. could be
#' represented by a column from a metadata file). Default = NULL.
#' @param covarDEA a character vector of covariate(s) to include in the analysis.
#' covarDEA can be defined only together with a batch parameter (which
#' represents one of the covariates). Default = NULL.
#' @param group a factor specifying group for each sample (e.g. could be
#' represented by a column from a metadata file)
#' @param cutoff.logFC A cutoff value for the logarithmic fold change applied
#' to each feature. Default = 1.
#' @param cutoff.FDR a cutoff value for the false discovery rate (corrected
#' p-value) applied to each feature. Default = 0.01.
#' @param prefix a character vector defining prefix of output file name.
#' @param block a factor specifying blocking variable. The block must be of
#' the same length as data and contain 2 or more options (e.g. could be
#' represented by a column from a metadata file). Default = NULL.
#' @export
#' @import limma
#' @import sva
#' @return a list of:
#' 1) re-formatted matrix of differential expression/abundance results from
#' limma.
#' 2) original results (a list) from limma
#' 3) a character vector with unique feature names
#' 4) a design matrix
#' 5) a contrast matrix
#' @examples {
#' ###run DEA
#' campp2_brca_1_DEA<-RunDEA(data=campp2_brca_1_normalized,
#' metadata=campp2_brca_1_meta,
#' group=campp2_brca_1_meta$subtype, prefix="test",
#' block=NULL, batch=campp2_brca_1_meta$age,
#' covarDEA = c("tumor_stage"), cutoff.logFC=1, cutoff.FDR=0.05)}


RunDEA <- function(data, metadata=NULL, group, batch=NULL, covarDEA=NULL, cutoff.logFC=1, cutoff.FDR=0.01, prefix, block=NULL) {

    ##test batch parameter
    if(!is.null(batch)){
        if (length(batch) != ncol(data)) {
            stop("Batch correction selected, but size of the batch metadata column does not match the number of samples!")
        }
    }

    ##test if both batch and covarDEA are defined
    if(is.null(batch) && !is.null(covarDEA)){
        stop("covarDEA cannot be defined without a batch parameter (representing one of the covariates.")
    }

    ####Test if the covarDEA value(s) is present in the metadata
    if(!is.null(covarDEA)){
        test<-NULL

        for(i in 1:length(covarDEA)){
            test[i] <- length(metadata[,colnames(metadata) %in% covarDEA[i]]) == nrow(metadata)
        }
        if(isFALSE(all(test))){  ##test if only TRUE is present
            stop("Please, check if covariates are defined properly for each sample - if a number of rows with covariates is equal to number of rows of the samples and/or if the column name in metadata table and the covarDEA agrees.")
        }
    }


    ###Create design matrix;f alternatively, Design.Matrix function from CAMPP2 can be used.
    if (is.null(batch) && is.null(covarDEA)) {
        design.matrix <- model.matrix(~0+group)
        out.name <- "_DEA"
    } else if (!is.null(batch) && is.null(covarDEA)) {
        design.matrix <- model.matrix(~0+group+batch)
        out.name <- "_databatch_DEA"
    } else if(!is.null(batch) && !is.null(covarDEA)) {
        covariates<-paste0(covarDEA, collapse = "+")
        design.matrix <- eval(parse(text=paste0("model.matrix(~0+group+batch","+",covariates,",data=metadata)")))  #This construct is able to process multiple covariates (multiple covariates not implemented in the Parse_arguments.R yet)
        out.name <- "_databatch_covars_DEA"
    }

    # Making group contrasts
    combinations <- data.frame(t(combn(paste0("group",levels(as.factor(group))), 2)))   ##we need this because "group" is implemented in the design matrix
    combinations$contr <- apply(combinations[,colnames(combinations)], 1, paste, collapse = "-")
    contrast.matrix <- makeContrasts(contrasts=combinations$contr,levels=as.character(colnames(design.matrix)))

    # Apply DEA to all comparisons
    res.DEA <- DEAFeatureApply(data, design.matrix, contrast.matrix, cutoff.logFC, cutoff.FDR, block)

    # Re-format res.DEA results and write results out as .txt file
    if (!is.null(res.DEA)) {
        DEA.out <- ExportDEA(res.DEA, paste0(prefix, out.name))   ###DEA.out is used as a main output
        rownames(DEA.out) <- NULL
        res.DEA.names <- unique(DEA.out$name)   ###could be used in LASSO
    } else {
        cat("No signficant DEA hits found. Check your cut-off for differential expression analysis, it may be that these are too stringent.")
    }


    return(list("DEA.out"=DEA.out,"res.DEA"=res.DEA,"res.DEA.names"=res.DEA.names, "DEA.design.matrix"=design.matrix, "DEA.contrast.matrix"=contrast.matrix)) ##the main output is DEA.out and res.DEA.names; keeping res.DEA this for possible troubleshooting in LASSO, WGCNA

}
