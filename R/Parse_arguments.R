#' @title Parse arguments
#' @description Parsing mandatory and optional arguments defined into the main function
#' @param data1 Gene count matrix (gene IDs as row names and sample IDs as columns). It's recommended to import gene counts using function "import_counts".
#' @param data2 Gene count matrix for a second dataset.
#' @param metadata1 Samples' metadata table should be imported using function "import_metadata". Metadata must include exactly the same samples as gene counts (data1) and samples must be sorted similarly.
#' @param metadata2 Metadata for a second dataset.
#' @param technology Technology used for the analysis of biological input. Current options are 'array', 'seq', 'ms' or 'other'. This argument is mandatory and depending on which option is chosen, data is transformed differently. If a second dataset is provided, the option should be specified for each dataset, provided as a character vector.
#' @param groups Argument defining groups of samples should be specified as a character vector. The first element specifying the name of the column in the metadata file containing sample IDs and the second element specifying the name of the column which contains the groups for the DE/DA analysis.
#' @param data.check Distributional checks of the input data is activated using logical argument (TRUE/FALSE). If activated, Cullen-Frey graphs will be made for 10 randomly selected variables to check data distributions. This argument is per default set to TRUE.
#' @param batches Specifies which metadata should be used for a batch correction (sequencing run/tissue/interstitial fluid/etc.). Argument takes a character vector of length 1 (one dataset) or 2 (two datasets), where the string(s) match a column name(s) in the metadata file(s). Default is NULL.
#' @param kmeans Argument specifies ("TRUE" or "FALSE") if a k-means clustering should be performed. Default is FALSE (do not run).
#' @param num.km.clusters a vector of manually defined number(s) of clusters. By default, the values(s) are calculated automatically (a default value is NULL).
#' @param plot.heatmap Argument for heatmap specified as either: "DE", "DA", "LASSO", "EN" or "Consensus". Defaults is FALSE (do not run).
#' @param correlation Argument for correlation analysis. String specify which features should be correlated, options are: "ALL", "DE", "DA", "LASSO", "EN" or "Consensus". For this type of analysis, 2 datasets must include the same samples, e.g. tumor1-normal vs tumor2-normal (3 samples from 1 patient needed). Default is FALSE (do not run).
#' @param survival (double check this when parsin survival function) Survival analysis may be performed on differentially expressed/abundant variables, variables from LASSO/EN regression or the consensus of these. Argument "survival" must be specified as either; "DE", "DA", "LASSO", "EN" or "Consensus". The full dataframe of variables may be used (if argument is set to ALL), HOWEVER this is not advisable unless the dataset is small with very few variables. At least, "survival", "outcome", "outcome.time" info must be included in the metadata file. The metadata file must contain at least four columns named; "ids" (sample identifiers), "age" (age in years at diagnosis, surgery or entry into trail), "outcome.time" (time until end of follow-up in weeks, months or years, censuring, death) and "outcome" (numeric 0 = censuring, 1=death). N.B. in case of (paired) normal samples the columns with survival information for these samples should contain "NA" values.
#' @param surv.plot Argument which specifies number of features to include per survival plot. Default is 50.
#' @param standardize (double check for sequencing) Data centering. This option may be set to "mean" or "median." If two datasets are provided, the standardize option should be specified for each dataset, provided as a character vector. If the argument standardize is not specified and "technology" = "array", then quantile normalization will be performed. Defaults is FALSE (do not run).
#' @param transform Data transformation type. Current options are "log2", "log10", "logit" and "voom". If two datasets are provided the parameter should be specified for each dataset, provided as a character vector. Defaults is FALSE (do not run).
#' @param prefix Prefix for the results' files and results folder. Defalt is "Results".
#' @param signif Cut-offs for log fold change (logFC) and corrected p-value (fdr), defining significant hits (proteins, genes, miRNAs or N-features). If argument is set, it must be a numeric vector, where the first element specifies the cut-off for logFC and the second element specifies the cut-off for corrected p-value (fdr).  In case of 2 datasets, vector must be of length 4. By default, cutoffs will be set to -1 > logFC > 1 and corrected p-values < 0.05.
#' @param plot.PCA This argument specifies ("TRUE" or "FALSE") if a preliminary PCAplot should be made for data overview. Default is FALSE (do not run).
#' @param show.PCA.labels a boolean value (TRUE or FALSE) specifying if elements (e.g. samples) should be labelled (for PCAPlot and runKmeans functions). Labeling is based on column names of the input data. Default value is FALSE.
#' @param covariates Covariates to include in the analysis. If multiple of these, they should be specified as a character vector. The first element in this list must be either TRUE or FALSE. If TRUE is specified then covariates will be included in both DE/DA analysis and Survival Analsysis. If FALSE is specified covariates will ONLY be used for Survival Analsysis. Names of covariates should match the desired columns in the metadata file. Default is NULL.
#' @param stratify This argument may be used if some of the categorical (NOT continous) covariates violate the cox proportional assumption. The workflow checks for proportional hazard and will retun the covariates that fail the PH test. You may then rerun the workflow with this argument followed by the names of the categorical covariates which failed and these will be stratified. Default is NULL.
#' @param block A vector or factor specifying a blocking variable for differential expression analysis. The block must be of same length as the data and contain 2 or more options. For 2 datasets, the block can be defined as a list of two vectors or factors.
#' @param colors Custom color pallet for PCA and heatmaps. Must be the same length as number of groups used for comparison (e.g. two groups = two colors) and must be defined as character vector. See R site for avalibe colors http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf. Default is NULL.
#' @param lasso Argument specifying parameters for LASSO or Elastic net regression. This argument may be set to 1 for LASSO or >0 & <1 for Elastic Net, but NOT to 0 exactly (Ridge Regression). Defaults is FALSE (do not run).
#' @param WGCNA Argument specifying parameter for Weighed Gene Co-expression Network Analysis. It takes a string, either "DA", "DE" or "ALL" specifying if all variables should be included in WGCNA or only differentially expressed / abundant variables. Defaults is FALSE (do not run).
#' @param cutoff.WGCNA Argument specifying the cutoff values for WGCNA. The argument takes a numuric vector of three values: (I) minimum modules size, (II) maximum % dissimilarity for merging of modules, and (III) % of top most interconnected genes (or other features) to return, from each modules identified in the Weighed Gene Co-expression Network Analysis. Default values are 10,25,25.
#' @param PPint Argument specifying that protein-protein interaction networks should be generated using the results of the differential expression analysis. This argument must be a character vector of length two. The first element in this list must be a string specifying the type of gene identifier in the gene counts file provided. Allowed identifiers are: "ensembl_peptide_id", "hgnc_symbol", "ensembl_gene_id", "ensembl_transcript_id", "uniprotswissprot". The second element is a string specifying version of the stringDB to use. Currently only version supported is: 11.0. Default is FALSE (do not run).
#' @param gene.miR.int Argument specifying that gene-miRNA interaction networks should be generated using the results of the differential expression analysis. This argument must be a character vector of length two. The first element in this list must be a string specifying the type of miRNA identifier in the gene counts data file. Allowed identifiers are: "mature_mirna_ids", "mature_mirna_accession". The second element must be a string specifying the miRNA-gene database to use, currently options are: "targetscan" (validated miRNAs), "mirtarbase" (predicted miRNAs), "tarscanbase" (validated + predicted miRNAs)". Default is FALSE (do not run).
#' @export
#' @return parsed arguments
#' @examples \dontrun{
#' ...
#' }

parseArguments <- function(data1, data2, metadata1, metadata2, groups, technology, batches, data.check, standardize, transform, plot.PCA, plot.heatmap, kmeans, num.km.clusters, signif, colors, block, prefix, correlation, lasso, WGCNA, cutoff.WGCNA, survival, covariates, stratify, surv.plot, PPint, gene.miR.int, show.PCA.labels){

    # For DE/DA analysis, survival analysis and correlation analysis
    DEA.allowed.type <- c("ALL","EN", "LASSO", "DA", "DE", "Consensus",FALSE)
    WGCNA.allowed.type <- c("DA", "DE", "ALL", FALSE)
    show.PCA.labels.allowed.type <- c(TRUE, FALSE)

    # For survival analysis, must be in metadata file:
    survival.metadata <- c("survival", "outcome", "outcome.time")

    # For miRNA-gene and protein-protein network analysis
    approved.gene.IDs <- c("ensembl_peptide_id", "hgnc_symbol","ensembl_gene_id","ensembl_transcript_id", "uniprotswissprot")
    approved.miR.IDs <- c("mature_mirna_ids", "mature_mirna_accession")
    gene.query <- c("stringdatabase")
    miR.query <- c("targetscan", "mirtarbase", "tarscanbase")


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

    # Groups
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

        # Groups
        group2 = as.factor(metadata2[ , groups[[4]]])
        if (length(group2) <= 1) {
            stop(paste0("No column in metadata2 file called ",groups[[4]]))
        }
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

    logFC=NULL
    slogFC=NULL
    FDR=NULL
    sFDR=NULL
    if (is.null(signif)){
        cat("\n- No cut-off for significant hits has been chosen. Cutoffs will be set to -1 > logFC > 1 and corrected p-value (fdr) < 0.05.")
        logFC <- 1
        FDR <- 0.05
        slogFC <- 1
        sFDR <- 0.05
    }
    else {
        signif <- as.numeric(signif)
        if (length(signif) == 2) {
            logFC <- signif[1]
            FDR <- signif[2]
        }
        else if (length(signif) == 4) {
            logFC <- signif[1]
            FDR <- signif[2]
            slogFC <- signif[3]
            sFDR <- signif[4]
        }
        else {
            stop("If argument signif is set, it must be a vector of length 2 OR 2*2 = 4 , if two datasets are used, (with quotes and parenthesis!) where the first element specifies the cut-off for logFC and the second element specifies the cut-off for corrected p-value (fdr) for each set. If signif is not specified defaults will be used. Cutoffs will be set to -1 > logFC > 1 and corrected p-value (fdr) < 0.05.")
        }
    }

    # Blocks

    block1=NULL
    block2=NULL
    if(!is.null(block)){
        if (class(block)=='list'){
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
        stop("Options for PCA labels are: TRUE or FALSE Please re-run pipeline with one of these!")
    }


    # Correlation

    if (correlation %in% c(DEA.allowed.type)) {
        corrby <- correlation
    } else {
        stop(paste0("Argument corr (correlation analysis) is specified. This argument takes a string denoting which set of variables to use for correlation analysis, options are: ", DEA.allowed.type,"."))
    }



    # WGCNA

    if (!WGCNA %in% WGCNA.allowed.type) {
        stop("WGCNA may be performed with either results of differential abundance / expression analysis (DA / DE) or with all variables (ALL). N.B It is not advisable to run WGCNA with all variables if n > 5000. This will make WGCNA slow and plots will be difficult to iterpret. If 'ALL' is chosen and n > 5000, NO plots will be generated, but module variable interconnectivity scores will still be computed.")
    }



    # CutoffWGCNA

    if (is.null(cutoff.WGCNA)){
        cutoff.WGCNA <- c(min(10, nrow(data1)/2), 25, 25)
    } else {
        cutoff.WGCNA <- as.numeric(cutoff.WGCNA)
    }



    # Survival Analysis

    if (!survival %in% c(DEA.allowed.type)) {
        stop("Options for survival analysis variable sets are: ALL,EN, LASSO, DA, DE, Consensus. Please re-run pipeline with one of these!")
    }



    # Covariates (DEA and survival)

    if (is.null(covariates)){
        covarDEA1 <- NULL
        covarDEA2 <- NULL
        covarS <- NULL
    } else {
        covarS <-  covariates[-1]
        if (covariates[1] == TRUE) {
            covarDEA1 <- covariates[-1]
            if (!is.null(data2)) {
                covarDEA2 <- covariates[-1]
            }
        } else if (covariates[1] == FALSE) {
            covarDEA1 <- NULL
            covarDEA2 <- NULL
        } else {
            stop("First argument in 'survival' must be TRUE or FALSE. If TRUE, covariates will be used for both DE analysis and survival analysis. If FALSE, covariates will be used only for survival analysis.")
        }
    }



    # Network Analysis
    PPI<-NULL
    if (!isFALSE(PPint)){
        PPI <- PPint
        if (!PPI[[1]] %in% approved.gene.IDs | !PPI[[2]] %in% gene.query | length(PPI) != 2) {
            stop(paste0("Argument x must be a vector of strings of legth two. First element in list must specify type of gene ID matching the type of IDs in expression file, approved options are: ", approved.gene.IDs, ".Second element must specify which PPI database to use, currently options are: ", gene.query, "."))
        }
    }

    GmiRI<-NULL
    if (!isFALSE(gene.miR.int)){
        GmiRI <- gene.miR.int
        if (!GmiRI[[1]] %in% approved.miR.IDs | !GmiRI[[2]] %in% miR.query | length(GmiRI) != 2) {
            stop(paste0("Argument x must be a vector of strings of legth two. First element in lit must specify type of miRNA ID matching the type of IDs in expression file, approved options are: ", approved.miR.IDs, ".Second element must specify which miRNA-gene database to use, currently options are: ", miR.query, "."))
        }
    }


    print("Printing defined/processed parameters:")
    cat(c("\n",
          paste0("technology: ",technology),"\n",
          paste0("groups: ",groups),"\n",
          #         paste0("group1: ",group1),"\n",
          #         paste0("group2: ",group2),"\n",
          #         paste0("ids: ",ids),"\n",
          paste0("batches: ",batches),"\n",
          #         paste0("batch1: ",batch1),"\n",
          #         paste0("batch2: ",batch2),"\n",
          paste0("databatch1: ",databatch1),"\n",
          paste0("databatch2: ",databatch2),"\n",
          paste0("standardize: ",standardize),"\n",
          paste0("transform: ",transform),"\n",
          paste0("data.check: ",data.check),"\n",
          paste0("plot.PCA: ",plot.PCA),"\n",
          paste0("kmeans: ",kmeans),"\n",
          paste0("num.km.clusters: ",num.km.clusters),"\n",
          paste0("signif: ",signif),"\n",
          paste0("logFC: ",logFC),"\n",
          paste0("FDR: ",FDR),"\n",
          paste0("slogFC: ",slogFC),"\n",
          paste0("sFDR: ",sFDR),"\n",
          paste0("blocks:", block),"\n",
          #         paste0("block1: ",block1),"\n",
          #         paste0("block2: ",block2),"\n",
          paste0("colors: ",colors),"\n",
          paste0("prefix: ",prefix),"\n",
          paste0("plot.heatmap: ",plot.heatmap),"\n",
          paste0("corrby: ",corrby),"\n",
          paste0("lasso: ",lasso),"\n",
          paste0("WGCNA: ",WGCNA),"\n",
          paste0("cutoff.WGCNA: ",cutoff.WGCNA),"\n",
          paste0("survival: ",survival),"\n",
          paste0("covarDEA1: ",covarDEA1),"\n",
          paste0("covarDEA2: ",covarDEA2),"\n",
          paste0("covarS: ",covarS),"\n",
          paste0("stratify: ",stratify),"\n",
          paste0("surv.plot: ",surv.plot),"\n",
          paste0("PPI: ",PPI),"\n",
          paste0("GmiRI: ",GmiRI),"\n",
          paste0("show.PCA.labels: ",show.PCA.labels),"\n"
    ))

    return(list("data1"=data1,"data2"=data2,"metadata1"=metadata1,"metadata2"=metadata2, "technology"=technology, "groups"=groups,"group1"=group1,"group2"=group2,"ids"=ids,"batches"=batches,"databatch1"=databatch1,"databatch2"=databatch2,"batch1"=batch1, "batch2"=batch2, "standardize"=standardize,"transform"=transform,"data.check"=data.check,"plot.PCA"=plot.PCA,"kmeans"=kmeans, "num.km.clusters"=num.km.clusters, "signif"=signif, "logFC"=logFC,"FDR"=FDR,"slogFC"=slogFC,"sFDR"=sFDR,"block"=block,"block1"=block1,"block2"=block2,"colors"=colors,"prefix"=prefix,"plot.heatmap"=plot.heatmap,"corrby"=corrby,"lasso"=lasso,"WGCNA"=WGCNA,"cutoff.WGCNA"=cutoff.WGCNA,"survival"=survival,"covarDEA1"=covarDEA1,"covarDEA2"=covarDEA2,"covarS"=covarS,"stratify"=stratify,"surv.plot"=surv.plot,"PPI"=PPI,"GmiRI"=GmiRI,"DEA.allowed.type"=DEA.allowed.type,"survival.metadata"=survival.metadata,"approved.gene.IDs"=approved.gene.IDs,"approved.miR.IDs"=approved.miR.IDs,"gene.query"=gene.query,"miR.query"=miR.query, "show.PCA.labels"=show.PCA.labels))

}

