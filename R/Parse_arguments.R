#' @title Parse arguments
#' @description Parsing mandatory and optional arguments defined into the main function
#' @param data Gene count matrix (genes as rows; samples as columns) should be imported using function "import_counts"; user can choose own way but has to keep gene IDs as row names and sample IDs as column names.
#' @param sdata Gene count matrix for a second dataset. Gene count matrix (genes as rows; samples as columns) should be imported using function "import_counts".
#' @param metadata Samples' metadata table, should be imported using function "import_metadata". Metadata must include exactly the same samples as gene counts (data) and samples must be sorted similarly.
#' @param smetadata Metadata for a second dataset, should be imported using function "import_metadata". Metadata to the secend dataset must include exactly the same samples as sdata and samples must be sorted similarly.
#' @param technology Technology used for the analysis of biological input. Current options are 'array', 'seq', 'ms' or 'other'. This argument is mandatory and depending on which option is chosen, data is transformed differently. If a second dataset is provided, the option should be specified for each dataset, provided as a character vector.
#' @param groups Argument defining groups of samples should be specified as a character vector. The first element specifying the name of the column in the metadata file containing sample IDs and the second element specifying the name of the column which contains the groups for the DE/DA analysis.
#' @param datacheck Distributional checks is defined using logical argument (TRUE/FALSE) which specifies whether Cullen-Frey graphs should be made for 10 randomly selected variables to check data distributions. This argument is per default set to TRUE.
#' @param batches Specifies which metadata should be used for a batch correction (sequencing run/tissue/interstitial fluid/etc.). Argument takes a character vector of length 1 (one dataset) or 2 (two datasets), where the string(s) match a column name(s) in the metadata file(s).
#' @param kmeans Argument for kmeans clustering. The parameter must be specified as a character vector matching the name of a column in the metadata file, denoting the labeling of points on the MDS plot(s). If a parameter is set to "TRUE" (no column name specified) no labels will be added to the plot. Works only for the first dataset (data).
#' @param plotheatmap Argument for heatmap specified as either: "DE", "DA", "LASSO", "EN" or "Consensus".
#' @param corr Argument for correlation analysis. String specify which features should be correlated, options are: "ALL", "DE", "DA", "LASSO", "EN" or "Consensus". For this type of analysis, 2 datasets must include the same samples, e.g. tumor1-normal vs tumor2-normal (3 samples from 1 patient needed).
#' @param survival Survival analysis may be performed on differentially expressed/abundant variables, variables from LASSO/EN regression or the consensus of these. Argument "survival" must be specified as either; "DE", "DA", "LASSO", "EN" or "Consensus". The full dataframe of variables may be used (if argument is set to ALL), HOWEVER this is not advisable unless the dataset is small with very few variables. Survival info must be included in the metadata file. The metadata file must contain at least four columns named; "ids" (sample identifiers), "age" (age in years at diagnosis, surgery or entry into trail), "outcome.time" (time until end of follow-up in weeks, months or years, censuring, death) and "outcome" (numeric 0 = censuring, 1=death). N.B. in case of (paired) normal samples the columns with survival information for these samples should contain "NA" values.
#' @param survplot Argument which specifies number of features to include per survival plot, e.g. many features requires splitting of the plot, default features per plot is 50.
#' @param standardize Data centering. This option may be set to "mean" or "median." If two datasets are provided, the standardize option should be specified for each dataset, provided as a character vector. If the argument standardize is not specified and "technology" = "array", then quantile normalization will be performed.
#' @param transform Data transformation type. Current options are "log2", "log10", "logit" and "voom". If two datasets are provided the parameter should be specified for each dataset, provided as a character vector. If argument is left out, no transformation of data will be performed.
#' @param prefix Prefix for the results' files and results folder.
#' @param signif Cut-offs for log fold change (logFC) and corrected p-value (fdr), defining significant hits (proteins, genes, miRNAs or N-features). If argument a signif is set, it must be a numeric vector, where the first element specifies the cut-off for logFC and the second element specifies the cut-off for corrected p-value (fdr).  In case of 2 datasets, vector must be of length 4. If omitted, cutoffs will be set to -1 > logFC > 1 and corrected p-values < 0.05.
#' @param plotmds This argument specifies ("TRUE" or "FALSE") if a preliminary MDSplot should be made for data overview. Works only for the first dataset.
#' @param covar Covariates to include in the analysis. If multiple of these, they should be specified as a character vector. The first element in this list must be either TRUE or FALSE. If TRUE is specified then covariates will be included in both DE/DA analysis and Survival Analsysis. If FALSE is specified covariates will ONLY be used for Survival Analsysis. Names of covariates should match the desired columns in the metadata file.
#' @param stratify This argument may be used if some of the categorical (NOT continous) covariates violate the cox proportional assumption. The pipline checks for proportional hazard and will retun the covariates that fail the PH test. You may then rerun the pipeline with this argument followed by the names of the categorical covariates which failed and these will be stratified.
#' @param colors Custom color pallet for MDS and heatmaps. Must be the same length as number of groups used for comparison (e.g. two groups = two colors) must be defined as character vector. See R site for avalibe colors http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf.
#' @param lasso Argument specifying parameters for LASSO or Elastic net regression. This argument may be set to 1 for LASSO or >0 & <1 for Elastic Net, but NOT to 0 exactly (Ridge Regression).
#' @param WGCNA Argument specifying parameter for Weighed Gene Co-expression Network Analysis. It takes a string, either "DA", "DE" or "ALL" specifying if all variables should be included in WGCNA or only differentially expressed / abundant variables.
#' @param cutoffWGCNA Argument specifying the cutoff values for WGCNA. The argument takes a numuric vector of three values: (I) minimum modules size, (II) maximum % dissimilarity for merging of modules, and (III) % of top most interconnected genes (or other features) to return, from each modules identified in the Weighed Gene Co-expression Network Analysis. Default values are 10,25,25.
#' @param PPint Argument specifying that protein-protein interaction networks should be generated using the results of the differential expression analysis. This argument must be a character vector of length two. The first element in this list must be a string specifying the type of gene identifier in the gene counts file provided. Allowed identifiers are: "ensembl_peptide_id", "hgnc_symbol", "ensembl_gene_id", "ensembl_transcript_id", "uniprotswissprot". The second element is a string specifying version of the stringDB to use. Currently only version supported is: 11.0.
#' @param GenemiRint Argument specifying that gene-miRNA interaction networks should be generated using the results of the differential expression analysis. This argument must be a character vector of length two. The first element in this list must be a string specifying the type of miRNA identifier in the gene counts data file. Allowed identifiers are: "mature_mirna_ids", "mature_mirna_accession". The second element must be a string specifying the miRNA-gene database to use, currently options are: "targetscan" (validated miRNAs), "mirtarbase" (predicted miRNAs), "tarscanbase" (validated + predicted miRNAs)".
#' @export
#' @return parsed arguments
#' @examples \dontrun{
#' ...
#' }

#Is there a better way how to create global variables?
parseArguments <- function(data, sdata, metadata, smetadata, groups, technology, batches, datacheck, standardize, transform, plotmds, plotheatmap, kmeans, signif, colors, prefix, corr, lasso, WGCNA, cutoffWGCNA, survival, covar, stratify, survplot, PPint, GenemiRint){

    # For DE/DA analysis, survival analysis and correlation analysis - user input must be in list:
    must.be.type <- c("ALL","EN", "LASSO", "DA", "DE", "Consensus")

    # For survival analysis, must be in metadata file:
    must.contain <- c("survival", "outcome", "outcome.time")

    # For miRNA-gene and protein-protein network analysis - user inputs must be in lists:
    approvedGeneIDs <- c("ensembl_peptide_id", "hgnc_symbol","ensembl_gene_id","ensembl_transcript_id", "uniprotswissprot")
    approvedmiRIDs <- c("mature_mirna_ids", "mature_mirna_accession")
    genequery <- c("stringdatabase")
    miRNAquery <- c("targetscan", "mirtarbase", "tarscanbase")

    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ### MANDATORY ARGUMENTS ###
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    if (is.null(technology)){
        stop("\n- Argument technology is missing. Technology specifies type of input data (and sdata, if included).\n")
    } else {
        technology <- technology
    }


    # IDs and Groups to check user input!
    if (is.null(groups)){
        stop("Argument groups should be specified as a vector of strings. The first element specifying the name of the column in the metadata file containing sample IDs and the second element specifying the name of the column which contains the groups for the DE/DA analysis.")
    } else {
        groups <- groups
        if (length(groups) < 2) {
            stop("Argument groups should be specified as a vector of strings. The first element specifying the name of the column in the metadata file containing sample IDs and the second element specifying the name of the column which contains the groups for the DE/DA analysis.")
        }
    }


    # IDs and Groups to contrast
    # IDs
    acall <- parse(text = paste0("metadata$", as.character(groups[[1]])))
    ids <- as.character(eval(acall))

    if (length(ids) <= 1) {
        stop(paste0("No column in metadata test ids file called ",groups[[1]]))
    } else {
        metadata$ids <- ids
    }


    # Match Data and Metadata
    metadata <- metadata[metadata$ids %in% gsub("^X","",colnames(data)),]  ###IMPORT OF XLSX adds "X" - MUST BE REMOVED (check.names cannot be deactivated)

    # Groups
    acall <- parse(text = paste0("metadata$", as.character(groups[[2]])))
    group <- as.factor(as.character(eval(acall)))

    if (length(group) <= 1) {
        stop(paste0("Check the column name for groups in your metadata ",groups[[2]], " If the name matches, check the samples names."))
    }


    sgroup=NULL
    if (length(groups) == 4) {

        # IDs
        acall <- parse(text = paste0("smetadata$", as.character(groups[[3]])))
        ids <- as.character(eval(acall))

        if (length(ids) <= 1) {
            stop(paste0("No column in metadata file called ",groups[[3]]))
        } else {
            smetadata$ids <- ids
        }

        # Match Data and Metadata
        smetadata <- smetadata[smetadata$ids %in% colnames(sdata),]

        # Groups
        acall <- parse(text = paste0("smetadata$", as.character(groups[[4]])))
        sgroup <- as.factor(as.character(eval(acall)))

        if (length(sgroup) <= 1) {
            stop(paste0("No column in metadata file called ",groups[[4]]))
        }
    }



    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ### ADDITIONAL ARGUMENTS ###
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



    # Batches
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    if (is.null(batches)){
        databatch <- FALSE
        sdatabatch <- FALSE
    } else {
        batches <- batches
        acall <- parse(text = paste0("metadata$", as.character(batches[[1]])))
        batch <- as.factor(as.character(eval(acall)))
        databatch <- TRUE
        if (length(batch) <= 1) {
            stop(paste0("No column in metadata file called ",as.character(batches[[1]])))
        }
        if (length(batches) > 1 & exists("smetadata")) {
            acall <- parse(text = paste0("smetadata$", as.character(batches[[2]])))
            sbatch <- as.factor(as.character(eval(acall)))
            sdatabatch <- TRUE
            if (length(sbatch) <= 1) {
                stop(paste0("No column in metadata file called ",as.character(batches[[2]])))
            }
        } else {
            sdatabatch <- FALSE
        }
    }



    # Standardize data
    if (is.null(standardize)){
        standardize <- c("none", "none")
    } else {
        standardize <- standardize
    }



    # Transform
    if (is.null(transform)){
        transform <- c("none", "none")
    } else {
        transform <- transform
    }


    # Check data distribution
    if (is.null(datacheck)){
        datacheck <- TRUE
    } else {
        datacheck <- datacheck
    }


    # MDS plot
    if (is.null(plotmds)){
        plotmds <- FALSE
    } else {
        plotmds <- plotmds
    }

    labels.kmeans=NULL
    # Kmeans
    if (is.null(kmeans)){
        kmeans <- FALSE
    } else {
        kmeans <- TRUE
        if (kmeans == TRUE) {
            labels.kmeans <- ""
        } else {
            file <- try(labels.kmeans <- as.character(eval(parse(text = paste0("metadata$", as.character(kmeans))))))
            if (class(file) == "try-error") {
                labels.kmeans <- ""
                rm(file)
            }
        }
    }



    # Significance
    if (is.null(signif)){
        cat("\n- No cut-off for significant hits has been chosen. Cutoffs will be set to -1 > logFC > 1 and corrected p-value (fdr) < 0.05.")
        logFC <- 1
        FDR <- 0.05
        slogFC <- 1
        sFDR <- 0.05
    } else {
        signif <- as.numeric(signif)
        if (length(signif) > 1) {
            logFC <- signif[1]
            FDR <- signif[2]
            if (length(signif) > 3) {
                slogFC <- signif[3]
                sFDR <- signif[4]
            }
        } else {
            stop("If argument signif is set, it must be a vector of length 2 OR 2*2 = 4 , if two datasets are used, (with quotes and parenthesis!) where the first element specifies the cut-off for logFC and the second element specifies the cut-off for corrected p-value (fdr) for each set. If signif is not specified defaults will be used. Cutoffs will be set to -1 > logFC > 1 and corrected p-value (fdr) < 0.05.")
        }
    }



    # Colors
    if (is.null(colors)){
        colors <- viridisLite::viridis(length(levels(group)), begin = 0.2, end = 0.8)
    } else {
        colors <- colors
    }


    # Prefix
    if (is.null(prefix)){
        prefix <- "Results"
    } else {
        prefix <- prefix
    }



    # Heatmap
    if (is.null(plotheatmap)){
        plotheatmap <- NULL
    } else {
        plotheatmap <- plotheatmap
    }



    # Correlation
    if (is.null(corr)){
        corrby <- NULL
    } else if (corr %in% must.be.type) {
        corrby <- corr
    } else {
        stop(paste0("Argument corr (correlation analysis) is specified. This argument takes a string denoting which set of variables to use for correlation analysis, options are: ", must.be.type,"."))
    }



    # LASSO
    if (is.null(lasso)){
        lasso <- NULL
    } else {
        lasso <- lasso
    }



    # WGCNA

    if (is.null(WGCNA)){
        WGCNA <- NULL
    } else if (WGCNA %in% c("DA", "DE", "ALL")) {
        WGCNA <- WGCNA
    } else {
        stop("WGCNA may be performed with either results of differential abundance / expression analysis (DA / DE) or with all variables (ALL). N.B It is not advisable to run WGCNA with all variables if n > 5000. This will make WGCNA slow and plots will be difficult to iterpret. If 'ALL' is chosen and n > 5000, NO plots will be generated, but module variable interconnectivity scores will still be computed.")
    }



    # CutoffWGCNA
    if (is.null(cutoffWGCNA)){
        cutoffWGCNA <- c(min(10, nrow(data)/2), 25, 25)
    } else {
        cutoffWGCNA <- as.numeric(cutoffWGCNA)
    }



    # Survival Analysis
    if (is.null(survival)){
        survival <- NULL
    } else {
        survival <- survival
        if ((!survival %in% must.be.type)) {
            stop("Options for survival analysis variable sets are: ALL,EN, LASSO, DA, DE, Consensus. Please re-run pipeline with one of these!")
        }
    }



    # Covariates (DEA and survival)
    if (is.null(covar)){
        covarD <- NULL
        scovarD <- NULL
        covarS <- NULL
    } else {
        covar <- covar
        covarS <-  covar[-1]
        if (covar[1] == TRUE) {
            covarD <- covar[-1]
            if (!is.null(sdata)) {
                scovarD <- covar[-1]
            }
        } else if (covar[1] == FALSE) {
            covarD <- NULL
            scovarD <- NULL
        } else {
            stop("First argument in 'survival' must be TRUE or FALSE. If TRUE, covariates will be used for both DE analysis and survival analysis. If FALSE, covariates will be used only for survival analysis.")
        }
    }



    # Stratify
    if (is.null(stratify)){
        stratify <- NULL
    } else {
        stratify <- stratify
    }



    # Survival Plot
    if (is.null(survplot)){
        survplot <- 50
    } else {
        survplot <- survplot
    }





    # Network Analysis

    if (is.null(PPint)){
        PPI <- NULL
    } else {
        PPI <- PPint
        if (!PPI[[1]] %in% approvedGeneIDs | !PPI[[2]] %in% genequery | length(PPI) != 2) {
            stop(paste0("Argument x must be a vector of strings of legth two. First element in list must specify type of gene ID matching the type of IDs in expression file, approved options are: ", approvedGeneIDs, ".Second element must specify which PPI database to use, currently options are: ", genequery, "."))
        }
    }


    if (is.null(GenemiRint)){
        GmiRI <- NULL
    } else {
        GmiRI <- GenemiRint
        if (!GmiRI[[1]] %in% approvedmiRIDs | !GmiRI[[2]] %in% miRNAquery | length(GmiRI) != 2) {
            stop(paste0("Argument x must be a vector of strings of legth two. First element in lit must specify type of miRNA ID matching the type of IDs in expression file, approved options are: ", approvedmiRIDs, ".Second element must specify which miRNA-gene database to use, currently options are: ", miRNAquery, "."))
        }
    }

    cat(c("\n",
         paste0("technology: ",technology),"\n",
         paste0("groups: ",groups),"\n",
         paste0("group: ",group[10,]),"\n",
         paste0("sgroup: ",sgroup),"\n",
         paste0("acall: ",acall),"\n",
#         paste0("ids: ",ids),"\n",
         paste0("batches: ",batches),"\n",
         paste0("databatch: ",databatch),"\n",
         paste0("sdatabatch: ",sdatabatch),"\n",
         paste0("standardize: ",standardize),"\n",
         paste0("transform: ",transform),"\n",
         paste0("datacheck: ",datacheck),"\n",
         paste0("plotmds: ",plotmds),"\n",
         paste0("kmeans: ",kmeans),"\n",
         paste0("labels.kmeans: ",labels.kmeans),"\n",
         paste0("signif: ",signif),"\n",
         paste0("logFC: ",logFC),"\n",
         paste0("FDR: ",FDR),"\n",
         paste0("slogFC: ",slogFC),"\n",
         paste0("sFDR: ",sFDR),"\n",
         paste0("colors: ",colors),"\n",
         paste0("prefix: ",prefix),"\n",
         paste0("plotheatmap: ",plotheatmap),"\n",
         paste0("corrby: ",corrby),"\n",
         paste0("lasso: ",lasso),"\n",
         paste0("WGCNA: ",WGCNA),"\n",
         paste0("cutoffWGCNA: ",cutoffWGCNA),"\n",
         paste0("survival: ",survival),"\n",
         paste0("covarD: ",covarD),"\n",
         paste0("scovarD: ",scovarD),"\n",
         paste0("covarS: ",covarS),"\n",
         paste0("stratify: ",stratify),"\n",
         paste0("survplot: ",survplot),"\n",
         paste0("PPI: ",PPI),"\n",
         paste0("GmiRI: ",GmiRI),"\n"
     ))
    my_list<-list("data"=data,"sdata"=sdata,"metadata"=metadata,"smetadata"=smetadata, "technology"=technology, "groups"=groups,"group"=group,"sgroup"=sgroup,"acall"=acall,"ids"=ids,"batches"=batches,"databatch"=databatch,"sdatabatch"=sdatabatch,"standardize"=standardize,"transform"=transform,"datacheck"=datacheck,"plotmds"=plotmds,"kmeans"=kmeans,"labels.kmeans"=labels.kmeans,"signif"=signif,"logFC"=logFC,"FDR"=FDR,"slogFC"=slogFC,"sFDR"=sFDR,"colors"=colors,"prefix"=prefix,"plotheatmap"=plotheatmap,"corrby"=corrby,"lasso"=lasso,"WGCNA"=WGCNA,"cutoffWGCNA"=cutoffWGCNA,"survival"=survival,"covarD"=covarD,"scovarD"=scovarD,"covarS"=covarS,"stratify"=stratify,"survplot"=survplot,"PPI"=PPI,"GmiRI"=GmiRI)
    return(my_list)
}
