#' @title Parse arguments
#' @description Parsing mandatory and optional arguments defined into the main function
#' @param data1 Gene count matrix (gene IDs as row names and sample IDs as columns). It's recommended to import gene counts using function "import_counts".
#' @param data2 Gene count matrix for a second dataset.
#' @param metadata1 Samples' metadata table should be imported using function "import_metadata". Metadata must include exactly the same samples as gene counts (data1) and samples must be sorted similarly.
#' @param metadata2 Metadata for a second dataset.
#' @param technology Technology used for the analysis of biological input. Current options are 'array', 'seq', 'ms' or 'other'. This argument is mandatory and depending on which option is chosen, data is transformed differently. If a second dataset is provided, the option should be specified for each dataset, provided as a character vector.
#' @param groups Argument defining groups of samples should be specified as a character vector. The first element specifying the name of the column in the metadata file containing sample IDs and the second element specifying the name of the column which contains the groups for the DE/DA analysis.
#' @param control.group A string vector defining control group name (e.g., "healthy" or "normal"). In case of subtype analysis (>2 groups), the output of the main wrapper will include comparisons between control group and each subtype.
#' @param data.check Distributional checks of the input data is activated using logical argument (TRUE/FALSE). If activated, Cullen-Frey graphs will be made for 10 randomly selected variables to check data distributions. This argument is per default set to TRUE.
#' @param plot.heatmap Argument defining which data will be used for the selection of the top x features to be plotted on the heatmap. Options are:"DEA", "LASSO", "EN", Ridge or "Consensus".
#' @param heatmap.size Argument specifying how many genes will be selected to be plotted on the heatmap if plot.heatmap is TRUE. The input must be specified as an even number.#' @param show.PCA.labels a boolean value (TRUE or FALSE) specifying if elements (e.g. samples) should be labelled (for PCAPlot and runKmeans functions). Labeling is based on column names of the input data. Default value is FALSE.
#' @param viridis.palette Argument specifying viridis color palette used for heatmaps. Default is "turbo".
#' @param batches Specifies which metadata should be used for a batch correction (sequencing run/tissue/interstitial fluid/etc.). Argument takes a character vector of length 1 (one data set) or 2 (two data sets), where the string(s) match a column name(s) in the metadata file(s). In case batch correction should be performed only in 1 out of 2 data sets, a data set without the batch correction (1st one in the example) should be define as "", e.g. batches(c("","column_name")). Default is NULL.
#' @param kmeans Argument specifies ("TRUE" or "FALSE") if a k-means clustering should be performed. Default is FALSE (do not run).
#' @param num.km.clusters either a vector of manually defined number(s) of clusters, or NULL. Each number in this vector represent a plausible number of clusters that the user would expected to be present in the data. If multiple values provided, the function will automatically perform K-means clustering using each of them as k (expected number of clusters) separately. If this argument is NULL (default), optimal numbers of clusters are calculated automatically based on the bayesian information criterion (mclust package), applied to sub sampled data (see documentation for the whole procedure).
#' @param standardize (double check for sequencing) Data centering. This option may be set to "mean" or "median." If two datasets are provided, the standardize option should be specified for each dataset, provided as a character vector. If the argument standardize is not specified and "technology" = "array", then quantile normalization will be performed. Defaults is FALSE (do not run).
#' @param transform Data transformation type. Current options are "log2", "log10", "logit" and "voom". If two datasets are provided the parameter should be specified for each dataset, provided as a character vector. Defaults is FALSE (do not run).
#' @param prefix Prefix for the results' files and results folder. Defalt is "Results".
#' @param signif Cut-offs for log fold change (logFC) and corrected p-value (fdr), defining significant hits (proteins, genes, miRNAs or N-features). If argument is set, it must be a numeric vector, where the first element specifies the cut-off for logFC and the second element specifies the cut-off for corrected p-value (fdr).  In case of 2 datasets, vector must be of length 4. By default, cutoffs will be set to -1 > logFC > 1 and corrected p-values < 0.05.
#' @param ensembl.version This argument specifies which ENSEMBL database to use when transforming ensemble IDs into HUGO IDs using biomaRt. The argument should be specified as a number.
#' @param plot.DEA This argument specifies ("TRUE" or "FALSE") whether visualizations should be made for the differential expression analysis.
#' @param plot.PCA This argument specifies ("TRUE" or "FALSE") if a preliminary PCAplot should be made for data overview. Default is FALSE (do not run).
#' @param show.PCA.labels a boolean value (TRUE or FALSE) specifying if elements (e.g. samples) should be labelled (for PCAPlot and runKmeans functions). Labeling is based on column names of the input data. Default value is FALSE.
#' @param covariates Covariates are specified as a character vector. Specified covariates will be included in both DEA analysis. Names of covariates should match the desired columns in the metadata file. Only one covariate for each dataset is allowed (multiple covariates are allowed when using RunDEA function out of this wrapper). Default is NULL.
#' @param block A vector or factor specifying a blocking variable for differential expression analysis. The block must be of same length as the data and contain 2 or more options. For 2 datasets, the block can be defined as a list of two vectors or factors.
#' @param colors Custom color pallet for MDS and heatmaps. Must be the same length as number of groups used for comparison (e.g. two groups = two colors) and must be defined as character vector. See R site for avalibe colors http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf. Default is NULL.
#' @param alpha.lasso a numeric vector specifying hyperparameter alpha for LASSO/Elastic network/Ridge regression. This value must be set to 0.0 < x < 1.0 for Elastic Net or to 1.0 for LASSO regression or to 0.0 for Ridge regression. Defaults is FALSE (do not run).
#' @param min.coef.lasso a numeric vector specifying a threshold for features' filtering (e.g. genes) based on the coefficients which are calculated during model fitting. Default value is > 0.
#' @param nfolds.lasso a numeric vector describing number of folds during Lambda estimation which is based on a cross-validation. Although nfolds can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. Smallest value allowable is nfolds=3. Default is 10.
#' @param num.trees.init an integer specifying number of trees to use for the first random forest in the random forest feature selection process. Default is NULL (not activated). Both num.trees.init and num.trees.iterat need to be > 0 to activate random forest.
#' @param num.trees.iterat an integer specifying number of trees to use for all additional random forests in the random forest feature selection process. Default is NULL (not activated). Both num.trees.init and num.trees.iterat need to be > 0 to activate random forest.
#' @param split.size an integer specifying the minimum number of samples that the groups must contain in order to carry out random forest classification and subsequent validation
#' @param test.train.ratio a floating point number between 0 and 1 representing the ratio of samples to keep as validation dataset. For example, a test.train.ratio = 0.25 splits 25 percent of the data into a validation dataset, meaning 75 percent of the data will be kept as the training dataset.
#' @export
#' @return parsed arguments
#' @examples \dontrun{
#' ...
#' }

parseArguments <- function(data1, data2, metadata1, metadata2, control.group, groups, technology, batches, data.check, standardize, transform, plot.PCA, plot.heatmap, plot.DEA, ensembl.version, heatmap.size, viridis.palette, kmeans, num.km.clusters, signif, colors, block, prefix, covariates, show.PCA.labels, alpha.lasso, min.coef.lasso, nfolds.lasso, num.trees.init, num.trees.iterat, split.size, test.train.ratio){


    # For DEA analysis, survival analysis and correlation analysis
    DEA.allowed.type <- c("ALL","EN", "LASSO", "Ridge", "DEA", "Consensus", FALSE)
    show.PCA.labels.allowed.type <- c(TRUE, FALSE)



        # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ### MANDATORY ARGUMENTS ###
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


    ### IDs and Groups to contrast
    # IDs
    ids = metadata1[ , groups[[1]]]

    if (length(ids) <= 1) {
        stop(paste0("No column in metadata1 ids file called ",groups[[1]]))
    } else {
        metadata1$ids <- ids
    }

    # groups
    group1 = as.factor(metadata1[ , groups[[2]]])

    if (length(group1) <= 1) {
        stop(paste0("Check the column name for groups in your metadata1 ",groups[[2]], " If the name matches, check the samples names."))
    }



    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ### OPTIONAL ARGUMENTS ###
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


    group2=NULL
    if (length(groups) == 4) {

        # IDs
        ids = metadata2[ , groups[[3]]]

        if (length(ids) <= 1) {
            stop(paste0("No column in metadata1 file called ",groups[[3]]))
        } else {
            metadata2$ids <- ids
        }


        # Match Data and Metadata
        metadata2 <- metadata2[metadata2$ids %in% colnames(data2),]

        # groups
        group2 = as.factor(metadata2[ , groups[[4]]])
        if (length(group2) <= 1) {
            stop(paste0("No column in metadata2 file called ",groups[[4]]))
        }

    }


    # control group
    if ((!isFALSE(plot.DEA) && is.null(control.group))){
        stop("For DEA fisualizations, control.group parameter must be defined.")
    }

    if(!is.null(control.group)){
        control.group1<-control.group[1]
        if (!is.null(control.group[2])){
            control.group2<-control.group[2]
        } else {
            control.group2<-NULL
        }
    } else {
        control.group1<-NULL
        control.group2<-NULL
    }


    # Batches

    batch1=NULL
    batch2=NULL
    if (is.null(batches)){
        databatch1 <- FALSE
        databatch2 <- FALSE
    } else {
        batch1 = as.factor(metadata1[ , batches[[1]]])
        databatch1 <- TRUE
        if (length(batch1) <= 1) {
            stop(paste0("No column in metadata1 file called ",as.character(batches[[1]])))
        }
        if (length(batches) > 1 & exists("metadata2")) {
            batch2 = as.factor(metadata2[ , batches[[2]]])
            databatch2 <- TRUE
            if (length(batch2) <= 1) {
                stop(paste0("No column in metadata2 file called ",as.character(batches[[2]])))
            }
        } else {
            databatch2 <- FALSE
        }
    }



    #Standardize

    if(isFALSE(standardize)){
        standardize<-c("none","none")
    }



    #Transform

    if(isFALSE(transform)){
        transform<-c("none","none")
    }



    # Significance

    cutoff.logFC1=NULL
    cutoff.logFC2=NULL
    cutoff.FDR1=NULL
    cutoff.FDR2=NULL
    if (is.null(signif)){
        cat("\n- No cut-off for significant hits has been chosen. Cutoffs will be set to -1 > logFC > 1 and corrected p-value (fdr) < 0.05.")
        cutoff.logFC1 <- 1
        cutoff.FDR1 <- 0.05
        cutoff.logFC2 <- 1
        cutoff.FDR2 <- 0.05
    }
    else {
        signif <- as.numeric(signif)
        if (length(signif) == 2) {
            cutoff.logFC1 <- signif[1]
            cutoff.FDR1 <- signif[2]
        }
        else if (length(signif) == 4) {
            cutoff.logFC1 <- as.numeric(signif[1])
            cutoff.FDR1 <- as.numeric(signif[2])
            cutoff.logFC2 <- as.numeric(signif[3])
            cutoff.FDR2 <- as.numeric(signif[4])
        }
        else {
            stop("If argument signif is set, it must be a vector of length 2 OR 2*2 = 4 , if two datasets are used, (with quotes and parenthesis!) where the first element specifies the cut-off for logFC and the second element specifies the cut-off for corrected p-value (fdr) for each set. If signif is not specified defaults will be used. Cutoffs will be set to -1 > logFC > 1 and corrected p-value (fdr) < 0.05.")
        }
    }

    # Blocks

    block1=NULL
    block2=NULL
    if(!is.null(block)){
        if(is(class(block), "list")){
            block1 <- block[[1]]
            block2 <- block[[2]]

            if ((length(block2)!=length(data2)) || (!length(unique(block2))>=2)){
                block2=NULL
                print("The given blocking parameter does not match the input data. Make sure that your input contains more than one block. Continuing analysis with standard value: block2=NULL")
            }

        } else {
            block1 <- block
        }
        if ((length(block1)!=length(data1)) || (!length(unique(block1))>=2)){
            block1=NULL
            print("The given blocking parameter does not match the input data. Make sure that your input contains more than one block. Continuing analysis with standard value: block1=NULL")
        }
    }

    # Colors

    if (is.null(colors)){
        colors <- viridisLite::viridis(length(unique(levels(c(group1,group2)))), begin = 0.2, end = 0.8)
    }


    #show.PCA.labels

    if (!show.PCA.labels %in% c(show.PCA.labels.allowed.type)) {
        stop("Allowed options for show.PCA.labes argument are: TRUE or FALSE. Please re-run pipeline with one of these!")
    }


    # Covariates (DEA)

    if (is.null(covariates)){
        covarDEA1 <- NULL
        covarDEA2 <- NULL
    } else {
        covarDEA1 <- covariates[1]
        if(!is.null(data2) && !is.null(covariates[2])){
            covarDEA2 <- covariates[2]
        }
    }

    print("Printing defined/processed parameters:")
    cat(c("\n",
          paste0("technology: ",technology),"\n",
          paste0("groups: ",groups),"\n",
          paste0("control.group1: ", control.group1), "\n",
          paste0("control.group2: ", control.group2), "\n",
          paste0("batches: ",batches),"\n",
          paste0("databatch1: ",databatch1),"\n",
          paste0("databatch2: ",databatch2),"\n",
          paste0("standardize: ",standardize),"\n",
          paste0("transform: ",transform),"\n",
          paste0("data.check: ",data.check),"\n",
          paste0("plot.PCA: ",plot.PCA),"\n",
          paste0("plot.DEA:",plot.DEA),"\n",
          paste0("kmeans: ",kmeans),"\n",
          paste0("num.km.clusters: ",num.km.clusters),"\n",
          paste0("cutoff.logFC1: ",cutoff.logFC1),"\n",
          paste0("cutoff.FDR1: ",cutoff.FDR1),"\n",
          paste0("cutoff.logFC2: ",cutoff.logFC2),"\n",
          paste0("cutoff.FDR2: ",cutoff.FDR2),"\n",
          paste0("blocks:", block),"\n",
          paste0("colors: ",colors),"\n",
          paste0("prefix: ",prefix),"\n",
          paste0("plot.heatmap: ",plot.heatmap),"\n",
          paste0("ensembl.version: ",ensembl.version),"\n",
          paste0("heatmap.size:",heatmap.size),"\n",
          paste0("viridis.palette:",viridis.palette), "\n",
          paste0("alpha.lasso: ",alpha.lasso),"\n",
          paste0("min.coef.lasso: ",min.coef.lasso),"\n",
          paste0("nfolds.lasso: ",nfolds.lasso),"\n",
          paste0("covarDEA1: ",covarDEA1),"\n",
          paste0("covarDEA2: ",covarDEA2),"\n",
          paste0("show.PCA.labels: ",show.PCA.labels),"\n",
          paste0("num.trees.init: ",num.trees.init),"\n",
          paste0("num.trees.iterat: ",num.trees.iterat),"\n",
          paste0("split.size: ",split.size),"\n",
          paste0("test.train.ratio: ",test.train.ratio),"\n"
    ))

    return(list("data1"=data1,"data2"=data2,"metadata1"=metadata1,"metadata2"=metadata2, "technology"=technology, "groups"=groups,
                "group1"=group1,"group2"=group2, 'control.group1' = control.group1,'control.group2' = control.group2, "ids"=ids,"batches"=batches,"databatch1"=databatch1,"databatch2"=databatch2,
                "batch1"=batch1, "batch2"=batch2, "standardize"=standardize,"transform"=transform,"data.check"=data.check,
                "plot.PCA"=plot.PCA,"kmeans"=kmeans, "num.km.clusters"=num.km.clusters, "cutoff.logFC1"=cutoff.logFC1,"cutoff.FDR1"=cutoff.FDR1,
                "cutoff.logFC2"=cutoff.logFC2,"cutoff.FDR2"=cutoff.FDR2,"block"=block,"block1"=block1,"block2"=block2,"colors"=colors,"prefix"=prefix,
                "plot.heatmap"=plot.heatmap,"covarDEA1"=covarDEA1,"covarDEA2"=covarDEA2,
                "DEA.allowed.type"=DEA.allowed.type,
                "show.PCA.labels"=show.PCA.labels,"heatmap.size"=heatmap.size,"viridis.palette"=viridis.palette,"ensembl.version"=ensembl.version, "plot.DEA"=plot.DEA,
                "alpha.lasso"=alpha.lasso,"min.coef.lasso"= min.coef.lasso,"nfolds.lasso"= nfolds.lasso, "num.trees.init"=num.trees.init, "num.trees.iterat"=num.trees.iterat, "split.size"=split.size, "test.train.ratio"=test.train.ratio))
}

