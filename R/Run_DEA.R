#' @title Run differential expression/abundance analysis
#' @descriptionA A function for running differential expression/abundance
#' analysis on a data matrix of features (e.g., genes) counts. Analysis is
#' based on voom transformation and limma.
#' @param data A raw gene count matrix from seq, array, ms or other technology
#' (with gene IDs as row names and sample IDs as columns). It's recommended to
#' import gene counts using function "import_counts".
#' @param metadata Samples' metadata table should be imported using function
#' "import_metadata". Metadata must include exactly the same samples as gene
#' counts (data) and samples must be sorted similarly. In case of this function,
#' metadata are used to extract covariates data for each sample. Default = NULL.
#' @param batch a factor specifying batch for each sample (e.g. could be
#' represented by a column from a metadata file). Default = NULL.
#' @param covarDEA a covariate(s) to include in the analysis provided
#' as a character vector. covarDEA can be defined only with a batch parameter
#' (which represents one of the covariates). Default = NULL.
#' @param group a factor specifying group for each sample (e.g. could be
#' represented by a column from a metadata file)
#' @param cutoff.logFC A cutoff value for the logarithmic fold change for each
#' feature. Default = 1.
#' @param cutoff.FDR a false discovery rate (corrected p-value) for each feature
#' @param prefix a character defining a prefix of output file name. Default =
#' 0.01.
#' @param block a factor specifying blocking variable. The block must be of
#' the same length as data and contain 2 or more options. Default = NULL.
#' @export
#' @import limma
#' @import sva
#' @seealso
#' @return a list of:
#' 1) re-formatted matrix of differential expression/abundance results
#' (containing information on gene regulation and significance of comparison)
#' 2) original results (a list) from limma
#' 3) a character vector with unique feature names
#' @examples \dontrun{
#' campp2_brca_1_DEA<-RunDEA(data=campp2_brca_1, metadata=campp2_brca_1_meta,
#' group=campp2_brca_1_meta$diagnosis, prefix="teset2",
#' block=campp2_brca_1_meta$subtype, batch=campp2_brca_1_meta$age,
#' covarDEA = c("tumor_stage"))


RunDEA <- function(data, metadata=NULL, group, batch=NULL, covarDEA=NULL, cutoff.logFC=1, cutoff.FDR=0.01, prefix, block=NULL) {

    ##test covarDEA parameter
    if(!is.null(covarDEA)){
        test<-NULL
        for(i in 1:length(covarDEA)){
            test[i] <- length(metadata[,colnames(metadata) %in% covarDEA[i]]) == nrow(metadata)
        }
        if(isFALSE(all(test))){  ##test if only TRUE is present
            stop("Please, check if covariates are defined properly for each sample - if a number of rows with covariates is equal to number of rows of the samples and/or if the column name in metadata table and the covarDEA agrees.")
        }
    }

    ##test batch parameter
    if(!is.null(batch)){
        if (length(batch) != ncol(data)) {
            stop("Batch correction selected, but size of the batch metadata column does not match the number of samples!")
        }
    }
    if(is.null(batch) && !is.null(covarDEA)){
        stop("covarDEA cannot be defined without a batch parameter (representing one of the covariates.")
    }

    ###Create design matrix; alternatively, Design.Matrix function can be used.
    if (is.null(batch) && is.null(covarDEA)) {
        design <- model.matrix(~0+group)
        out.name <- "_DEA"
    } else if (!is.null(batch) && is.null(covarDEA)) {
        design <- model.matrix(~0+group+batch)
        out.name <- "_databatch_DEA"
    } else if(!is.null(batch) && !is.null(covarDEA)) {
        names<-paste0(covarDEA, collapse = "+")
        design <- eval(parse(text=paste0("model.matrix(~0+group+batch","+",names,",data=metadata)")))
        out.name <- "_databatch_covars_DEA"
    }

    # Making group contrasts
    combinations <- data.frame(t(combn(paste0("group",levels(as.factor(group))), 2)))   ##we need this because "group" is implemented in the design matrix
    combinations$contr <- apply(combinations[,colnames(combinations)], 1, paste, collapse = "-")
    contrast.matrix <- makeContrasts(contrasts=combinations$contr,levels=as.character(colnames(design)))

    # Apply DEA to all comparisons
    res.DEA <- DEAFeatureApply(data, design, contrast.matrix, cutoff.logFC, cutoff.FDR, block)
    print("res.DEA")
    print(res.DEA)

    # Re-format res.DEA results and write results out as .txt file
    if (!is.null(res.DEA)) {
        DEA.out <- ExportDEA(res.DEA, paste0(prefix, out.name))   ###DEA.out is used as a main output
        rownames(DEA.out) <- NULL
        res.DEA.names <- unique(DEA.out$name)   ###could be used in LASSO
    } else {
        cat("No signficant DEA hits found. Check your cut-off for differential expression analysis, it may be that these are too stringent.")
    }

    return(list("DEA.out"=DEA.out,"res.DEA"=res.DEA,"res.DEA.names"=res.DEA.names)) ##keeping res.DEA this for possible troubleshooting in LASSO, WGCNA

}
