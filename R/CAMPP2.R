#' @title CAMPP2 pipeline
#' @description CAMPP2 is a tool for quantitative data analysis
#' @param data1 Gene count matrix (gene IDs as row names and sample IDs as columns). It's recommended to import gene counts using function "import_counts".
#' @param data2 Gene count matrix for a second dataset.
#' @param metadata1 Samples' metadata table should be imported using function "import_metadata". Metadata must include exactly the same samples as gene counts (data1) and samples must be sorted similarly.
#' @param metadata2 Metadata for a second dataset.
#' @param technology Technology used for the analysis of biological input. Current options are 'array', 'seq', 'ms' or 'other'. This argument is mandatory and depending on which option is chosen, data is transformed differently. If a second dataset is provided, the option should be specified for each dataset, provided as a character vector.
#' @param groups Argument defining groups of samples should be specified as a character vector. The first element specifying the name of the column in the metadata file containing sample IDs and the second element specifying the name of the column which contains the groups for the DEA analysis.
#' @param data.check Distributional checks of the input data is activated using logical argument (TRUE/FALSE). If activated, Cullen-Frey graphs will be made for 10 randomly selected variables to check data distributions. This argument is per default set to TRUE.
#' @param batches Specifies which metadata should be used for a batch correction (sequencing run/tissue/interstitial fluid/etc.). Argument takes a character vector of length 1 (one dataset) or 2 (two datasets), where the string(s) match a column name(s) in the metadata file(s). Default is NULL.
#' @param kmeans Argument for kmeans clustering. The parameter must be specified as a character vector matching the name of a column in the metadata file, denoting the labeling of points on the MDS plot(s). If a parameter is set to "TRUE" (no column name specified) no labels will be added to the plot. Works only for the first dataset (data1). Default is FALSE (do not run).
#' @param plot.heatmap Argument for heatmap specified as either: "DEA", "LASSO", "EN" or "Consensus". Defaults is FALSE (do not run).
#' @param heatmap.size Argument specifying how many genes will be plotted in the heatmap if plot.heatmap is TRUE. The input must be specified as an even number.
#' @param correlation Argument for correlation analysis. String specify which features should be correlated, options are: "ALL", "DEA", "LASSO", "EN" or "Consensus". For this type of analysis, 2 datasets must include the same samples, e.g. tumor1-normal vs tumor2-normal (3 samples from 1 patient needed). Default is FALSE (do not run).
#' @param survival (double check this when parsin survival function) Survival analysis may be performed on differentially expressed/abundant variables, variables from LASSO/EN regression or the consensus of these. Argument "survival" must be specified as either; "DEA", "LASSO", "EN" or "Consensus". The full dataframe of variables may be used (if argument is set to ALL), HOWEVER this is not advisable unless the dataset is small with very few variables. At least, "survival", "outcome", "outcome.time" info must be included in the metadata file. The metadata file must contain at least four columns named; "ids" (sample identifiers), "age" (age in years at diagnosis, surgery or entry into trail), "outcome.time" (time until end of follow-up in weeks, months or years, censuring, death) and "outcome" (numeric 0 = censuring, 1=death). N.B. in case of (paired) normal samples the columns with survival information for these samples should contain "NA" values.
#' @param surv.plot Argument which specifies number of features to include per survival plot. Default is 50.
#' @param standardize (double check for sequencing) Data centering. This option may be set to "mean" or "median." If two datasets are provided, the standardize option should be specified for each dataset, provided as a character vector. If the argument standardize is not specified and "technology" = "array", then quantile normalization will be performed. Defaults is FALSE (do not run).
#' @param transform Data transformation type. Current options are "log2", "log10", "logit" and "voom". If two datasets are provided the parameter should be specified for each dataset, provided as a character vector. Defaults is FALSE (do not run).
#' @param prefix Prefix for the results' files and results folder. Defalt is "Results".
#' @param signif Cut-offs for log fold change (logFC) and corrected p-value (fdr), defining significant hits (proteins, genes, miRNAs or N-features). If argument is set, it must be a numeric vector, where the first element specifies the cut-off for logFC and the second element specifies the cut-off for corrected p-value (fdr).  In case of 2 datasets, vector must be of length 4. By default, cutoffs will be set to -1 > logFC > 1 and corrected p-values < 0.05.
#' @param plot.mds This argument specifies ("TRUE" or "FALSE") if a preliminary MDSplot should be made for data overview. Works only for the first dataset. Default is FALSE (do not run).
#' @param MDS.labels Argument specifying whether to include labels on the MDS plot if plot.mds is TRUE. The argument should be defined with TRUE/FALSE.
#' @param plot.DEA This argument specifies ("TRUE" or "FALSE") whether visualizations should be made for the differential expression analysis. Visualizations include a Volcano plot for the groups of samples and an Upset plot and Venn diagram for sample subtypes in the input data if subtypes are present.
#' @param ensembl.version This arugment specifies which ensembl database to use when transforming ensemble IDs into HUGO IDs using biomaRt. The argument should be specified as a number.
#' @param covariates Covariates to include in the analysis. If multiple of these, they should be specified as a character vector. The first element in this list must be either TRUE or FALSE. If TRUE is specified then covariates will be included in both DE/DA analysis and Survival Analsysis. If FALSE is specified covariates will ONLY be used for Survival Analsysis. Names of covariates should match the desired columns in the metadata file. Default is NULL.
#' @param stratify This argument may be used if some of the categorical (NOT continous) covariates violate the cox proportional assumption. The workflow checks for proportional hazard and will retun the covariates that fail the PH test. You may then rerun the workflow with this argument followed by the names of the categorical covariates which failed and these will be stratified. Default is NULL.
#' @param colors Custom color pallet for MDS and heatmaps. Must be the same length as number of groups used for comparison (e.g. two groups = two colors) and must be defined as character vector. See R site for avalibe colors http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf. Default is NULL.
#' @param block A vector or factor specifying a blocking variable for differential expression analysis. The block must be of same length as the belonging dataset and contain 2 or more options. For 2 datasets the block can be defined as a list of factors or vectors.
#' @param lasso Argument specifying parameters for LASSO or Elastic net regression. This argument may be set to 1 for LASSO or >0 & <1 for Elastic Net, but NOT to 0 exactly (Ridge Regression). Defaults is FALSE (do not run).
#' @param WGCNA Argument specifying parameter for Weighed Gene Co-expression Network Analysis. It takes a string, either "DEA" or "ALL" specifying if all variables should be included in WGCNA or only differentially expressed / abundant variables. Defaults is FALSE (do not run).
#' @param cutoff.WGCNA Argument specifying the cutoff values for WGCNA. The argument takes a numuric vector of three values: (I) minimum modules size, (II) maximum percentage of dissimilarity for merging of modules, and (III) percentage of top most interconnected genes (or other features) to return, from each modules identified in the Weighed Gene Co-expression Network Analysis. Default values are 10,25,25.
#' @param PPint Argument specifying that protein-protein interaction networks should be generated using the results of the differential expression analysis. This argument must be a character vector of length two. The first element in this list must be a string specifying the type of gene identifier in the gene counts file provided. Allowed identifiers are: "ensembl_peptide_id", "hgnc_symbol", "ensembl_gene_id", "ensembl_transcript_id", "uniprotswissprot". The second element is a string specifying version of the stringDB to use. Currently only version supported is: 11.0. Default is FALSE (do not run).
#' @param gene.miR.int Argument specifying that gene-miRNA interaction networks should be generated using the results of the differential expression analysis. This argument must be a character vector of length two. The first element in this list must be a string specifying the type of miRNA identifier in the gene counts data file. Allowed identifiers are: "mature_mirna_ids", "mature_mirna_accession". The second element must be a string specifying the miRNA-gene database to use, currently options are: "targetscan" (validated miRNAs), "mirtarbase" (predicted miRNAs), "tarscanbase" (validated + predicted miRNAs)". Default is FALSE (do not run).
#' @import zeallot
#' @export
#' @seealso
#' @return CAMPP2 results
#' @examples \dontrun{
#' runCampp2(batches=c("tumor_stage","tumor_stage"),prefix="test_CAMPP2", data1=dataset1, data2=dataset2, metadata1=metadata1,metadata2=metadata2, groups=c("IDs", "diagnosis","IDs", "diagnosis"), technology=c("seq","seq"))
#' }


runCampp2 <- function (data1, metadata1, data2=NULL, metadata2=NULL, technology, groups, batches=NULL, data.check=TRUE, standardize=FALSE, transform=FALSE, plot.mds=FALSE, MDS.labels=TRUE, plot.DEA=FALSE, plot.heatmap=FALSE, heatmap.size=40, ensembl.version=104, kmeans=FALSE, signif=NULL, block=NULL, colors=NULL, prefix="Results", correlation=FALSE, lasso=FALSE, WGCNA=FALSE, cutoff.WGCNA=NULL, survival=FALSE, covariates=NULL, stratify=NULL, surv.plot=50, PPint=FALSE, gene.miR.int=FALSE){


    ###parse input arguments and assign updated values
    c(data1,data2,metadata1,metadata2,technology,groups,
      group1,group2,ids,batches,databatch1,databatch2,
      batch1,batch2,standardize,transform,data.check,
      plot.mds,MDS.labels,kmeans,labels.kmeans,signif,logFC1,FDR1,
      logFC2,FDR2,block,block1,block2, colors,prefix,plot.DEA, plot.heatmap,
      heatmap.size,ensembl.version,corrby,lasso,WGCNA,cutoff.WGCNA,
      survival,covarDEA1,covarDEA2,covarS,stratify,surv.plot,PPI,GmiRI,
      DEA.allowed.type,survival.metadata,approved.gene.IDs,provedmiRIDs,gene.query,miR.query) %<-% parseArguments(data1=data1, metadata1=metadata1, data2=data2,
                                                                                                    metadata2=metadata2,groups=groups, technology=technology,prefix=prefix,
                                                                                                    batches=batches,data.check=data.check,standardize=standardize, transform=transform,
                                                                                                    plot.mds=plot.mds,MDS.labels=MDS.labels, plot.DEA=plot.DEA, plot.heatmap=plot.heatmap,
                                                                                                    heatmap.size=heatmap.size,ensembl.version=ensembl.version,kmeans=kmeans,signif=signif,
												                                                    block=block, colors=colors,correlation=correlation,lasso=lasso,WGCNA=WGCNA,
												                                                    cutoff.WGCNA=cutoff.WGCNA,survival=survival,covariates=covariates,
												                                                    stratify=stratify,surv.plot=surv.plot,PPint=PPint, gene.miR.int=gene.miR.int)


    # Create directory Results

    dir.create(prefix)
    setwd(paste0(prefix, "/"))


    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #                                                                         ## MISSING VALUE IMPUTATIONS ###
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    print("RUNNING MISSING VALUE IMPUTATIONS")

    print("Running missing values imputation on data1")
    data1<-ReplaceNAs(data1)
    print("Missing values imputation on data1 has finished")

    if (!is.null(data2)){
        print("Running missing values imputation on data2")
        data2<-ReplaceNAs(data2)
        print("Missing values imputation on data2 has finished")
    }

    print("MISSING VALUE IMPUTATIONS FINISHED")

    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #                                                                         ## Fix zeros and check for negative values. ###
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    print("CAMPP2 automatically detects negative values and fix zeros in your data")
    print("RUNNING FIXING OF NEGATIVE AND ZERO VALUES")
    print("Detecting negative values and fixing zeros in data1")

    data1.original <- data1
    data1 %<-% FixZeros(data1,group1)

    data2.original=NULL
    if (!is.null(data2)){
        data2.original <- data2
        print("Detecting negative values and replacing zeros in data2")
        data2 %<-% FixZeros(data2,group2)
    }

    print("FIXING OF NEGATIVE AND ZERO VALUES FINISHED")

    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #                                                                         ## Normalization, Filtering and Transformation ###
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    print("CAMPP2 automatically performs data normalization and transformation depending on technology from which data are derived.")
    print("PROCESSING NORMALIZATION AND TRANSFORMATION")
    print("Normalizing and tranforming data1")

    if(!is.null(data1.original) && technology[1] == "seq"){    ###Only sequencing data could be normalized with 0-values
        print("Original data1 including 0-values are being normalized")
        data1 <- NormalizeData(data1.original, group1, transform[1], standardize[1], technology[1])
    } else {
        print("data1 without 0-values are being normalized")
        data1 <- NormalizeData(data1, group1, transform[1], standardize[1], technology[1])
    }

    if(!is.null(data2) || !is.null(data2.original)){
        print("Normalizing and tranforming data2")
    if (length(technology) != 2 || length(technology) == 1) {
        stop("\nTwo datasets are defined in the analysis, BUT argument technology has defined only one. Technology must be defined as string vector of length two.\n")
    } else if (length(transform) != 2 || length(transform) == 1) {
        stop("\nTwo datasets are defined in the analysis, BUT argument transform has defined only one. Transform must be defined as string vector of length two.\n")
    }

    if(!is.null(data2.original) && technology[2] == "seq"){       ###Only sequencing data could be normalized with 0-values
        print("Original data2 including 0-values are being normalized")
        data2 <- NormalizeData(data2.original, group2, transform[2], standardize[2], technology[2])
    } else {
        print("data2 without 0-values are being normalized")
        data2 <- NormalizeData(data2, group2, transform[2], standardize[2], technology[2])
        }
    }


    print("PROCESSING NORMALIZATION AND TRANSFORMATION FINISHED")


    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ### BATCH CORRECTION ###
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    print("RUNNING BATCH CORRECTION")


    if (databatch1==TRUE){
        print("Run batch correction on the 1st dataset")
        data1.batch %<-% BatchCorrect(data1,batch1,group1,technology[1])
        print("Batch correction of the first dataset finished")
    }else{
        print("Batch correction wasn't selected")
    }

    if (databatch2==TRUE){
        print("Run batch correction on the 2nd dataset")
        data2.batch %<-% BatchCorrect(data2,batch2,group2,technology[2])
        print("Batch correction of the second dataset finished")
    }else{
        print("Batch correction wasn't selected")
    }

    print("BATCH CORRECTION PART FINISHED")

    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ### Distributional Checks ###
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    print("PROCESSING DISTRIBUTIONAL CHECK")

    if (data.check == TRUE) {
        if (databatch1 == TRUE) {
            subset.data1 <- data1.batch[sample(nrow(data1.batch), 10),]
        } else {
            if (technology[1] == "seq") {
                subset.data1 <- data1$E[sample(nrow(data1$E), 10),]
            } else {
                subset.data1 <- data1[sample(nrow(data1), 10),]
            }
        }

        list.of.dists <- FitDistributions(subset.data1)

        dir.create("DataChecks")
        setwd("DataChecks/")
        PlotDistributions(subset.data1, list.of.dists)
        setwd("..")

        rm(subset.data1, list.of.dists)
    }

    if (data.check == TRUE & !is.null(data2)) {
        if (databatch2 == TRUE) {
            subset.data1 <- data2.batch[sample(nrow(data2.batch), 10),]
        } else {

            if (technology[2] == "seq") {
                subset.data1 <- data2$E[sample(nrow(data2$E), 10),]
            } else {
                subset.data1 <- data2[sample(nrow(data2), 10),]
            }
        }

        list.of.dists <- FitDistributions(subset.data1)

        dir.create("SecondDataChecks")
        setwd("SecondDataChecks/")
        PlotDistributions(subset.data1, list.of.dists)
        setwd("..")

        rm(subset.data1, list.of.dists)
    }

    print("DISTRIBUTIONAL CHECK PART FINISHED")

    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ### PRELIMINARY MDS PLOT ###
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    print("PROCESSING PRELIMINARY MDS PLOT")

    MDScolors <- gsub(pattern = "FF", replacement = "", x = colors)


    if (plot.mds == TRUE && databatch1 == TRUE){
        mdsplot <- MDSPlot(data1.batch, group1, ids, MDScolors, MDS.labels)
        ggsave(paste0(prefix, "_MDSplot_batchcorr.pdf"), plot = mdsplot, dpi = 300, width = 8, height = 8)

    } else if (plot.mds == TRUE && databatch1 == FALSE){
        if (technology[1] == "seq") {
            mdsplot <- MDSPlot(data.frame(data1$E), group1, ids, MDScolors, MDS.labels)
        } else {
            mdsplot <- MDSPlot(data1, group1, ids, MDScolors, MDS.labels)
        }
        ggsave(paste0(prefix, "_MDSplot.pdf"), plot = mdsplot, dpi = 300, width = 8, height = 8)
        rm(mdsplot)

    } else {
        cat("\n- No preliminary plot requested.\n")
    }

    print("PRELIMINARY MDS PLOT PART FINISHED")


    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ### Kmeans ###
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    print("PROCESSING KMEANS")

    if (kmeans == TRUE) {

        if (technology[1] == "seq") {
            k.data1 <- data1$E
        } else {
            k.data1 <- data1
        }

        # Number of sample sets to generate
        nsets <- 1:ceiling(nrow(k.data1)/1000)
        if(length(nsets) > 10) {
            nsets <- 1:10
        }

        # Number of variables in each sample set
        setsize <- nrow(k.data1)
        if (setsize > 2000) {
            setsize <- 2000
        }


        # Number of kmeans to try
        if(ncol(k.data1) <= 100) {
            nks <- 2:6
        } else if (ncol(k.data1) > 100 && ncol(k.data1) <= 500) {
            nks <- 2:11
        } else {
            nks <- 2:16
        }

        cat(paste0("Based on size of dataset, ", length(nsets), " sample sets will be generated of size ", setsize, " and ", length(nks), " clusters will be tested - N.B This may take up to 10 min!\nRunning......"))

        dir.create("KmeansResults")
        setwd("KmeansResults/")


        list.of.dfs <- list()

        if (databatch1 == TRUE) {
            for (idx in 1:length(nsets)) {
                df <- t(data1.batch[sample(nrow(data1.batch), setsize), ])
                list.of.dfs[[idx]] <- df
            }
            Kmeans.list <- lapply(list.of.dfs, function(x) EstimateKmeans(x, nsets))
            Kmeans.Out <- PlotKmeans(data1.batch, Kmeans.list, nks, labels.kmeans, prefix)
        } else {
            for (idx in 1:length(nsets)) {
                df <- t(k.data1[sample(nrow(k.data1), setsize), ])
                list.of.dfs[[idx]] <- df
            }
            Kmeans.list <- lapply(list.of.dfs, function(x) EstimateKmeans(x, nsets))
            Kmeans.Out <- PlotKmeans(k.data1, Kmeans.list, nks, labels.kmeans, prefix)
        }


        out <- cbind(metadata1, Kmeans.Out)

        write.table(out, paste0(prefix,"_Metadata_Kmeans.txt"), sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

        setwd("..")

        rm(k.data1, Kmeans.list, Kmeans.Out, out)
    }


    print("KMEANS PART FINISHED")


    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # CLEAN UP AND SOURCE FUNCTIONS FROM THE FUNCTIONS SCRIPT - PART 2
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    #rm(SplitList, ReadMyFile, ReplaceNAs, ReplaceZero, NormalizeData, FitDistributions, PlotDistributions, MDSPlot, EstimateKmeans)
    gc(full = TRUE)



    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ### DIFFERENTIAL EXPRESSION ANALYSIS ###
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    print("PROCESSING DIFFERENTIAL EXPRESSION")

    # Differential Expression Analysis with Limma
    dir.create("Results_DEA")
    setwd("Results_DEA/")

    #First dataset
    DEARes1 <- RunDEA(data1, technology[1], batch1, covarDEA1, group1, logFC1, FDR1, paste0(prefix,"1"), block1)

    DEA1.out<-DEARes1$DEA.out

    if (length(unique(DEA1.out$comparison)) > 1){
        DEA1.out<-subset(DEA1.out, grepl("healthy", comparison, fixed = TRUE))
    }

    res.DEA1<-DEARes1$res.DEA
    res.DEA1.names<-DEARes1$res.DEA.names

    #Second dataset
    if(!is.null(data2) & !is.null(metadata2)) {
        DEARes2 <- RunDEA(data2, technology[2], batch2, covarDEA2, group2, logFC2, FDR2, paste0(prefix,"2"), block2)

        DEA2.out<-DEARes2$DEA.out

        if (length(unique(DEA2.out$comparison)) > 1){
            DEA2.out<-subset(DEA2.out, grepl("healthy", comparison, fixed = TRUE))
        }

        res.DEA2<-DEARes2$res.DEA
        res.DEA2.names<-DEARes2$res.DEA.names
    }

    setwd("..")

    print("DIFFERENTIAL EXPRESSION PART DONE")

    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ### VISUALISATIONS FOR DIFFERENTIAL EXPRESSION ANALYSIS ###
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    setwd("Results_DEA/")

    if (!isFALSE(plot.DEA)){
        print("PROCESSING VISUALISATIONS FOR DIFFERENTIAL EXPRESSION ANALYSIS")

        #First dataset
        DEA.visuals <- AddGeneName(DEA1.out,ensembl.version)
        MakeVolcano(DEA.visuals, prefix,logFC1,FDR1)

        #Subtype analysis
        if (length(unique(DEA.visuals$comparison)) > 1){
            subtypes <- split.data.frame(DEA.visuals,DEA.visuals$comparison)

            sublist=list()
            for (genes in subtypes){
                sublist <- append(sublist,list(genes$Ensembl_ID))
            }
            MakeUpset(prefix,sublist,names(subtypes))
            MakeVennDiagram(prefix,sublist,names(subtypes))

            #Subtypes and expression
            subtypes_updown<-split.data.frame(DEA.visuals,list(DEA.visuals$comparison,DEA.visuals$dir))
            sublist_updown=list()
            for (genes in subtypes_updown){
                sublist_updown <- append(sublist_updown,list(genes$Ensembl_ID))
            }

            MakeUpset(paste0(prefix,".updown"),sublist_updown,names(subtypes_updown))
        }

        #Second dataset
        if (!is.null(data2) & !is.null(metadata2)){
            DEA2.visuals <- AddGeneName(DEA2.out,ensembl.version)
            MakeVolcano(DEA2.visuals, paste0(prefix,"_second"),logFC2,FDR2)

            #Subtype analysis
            if (length(unique(DEA2.visuals$comparison)) > 1){
                subtypes <- split.data.frame(DEA2.visuals,DEA2.visuals$comparison)

                sublist=list()
                for (genes in subtypes){
                    sublist <- append(sublist,list(genes$Ensembl_ID))
                }
                MakeUpset(paste0(prefix,"_second"),sublist,names(subtypes))
                MakeVennDiagram(paste0(prefix,"_second"),sublist,names(subtypes))

                #Subtypes and expression
                subtypes_updown<-split.data.frame(DEA2.visuals,list(DEA2.visuals$comparison,DEA2.visuals$dir))
                sublist=list()
                for (genes in subtypes_updown){
                    sublist <- append(sublist,list(genes$Ensembl_ID))
                }
                MakeUpset(paste0(prefix,".updown_second"),sublist,names(subtypes_updown))
            }
        }
    print("VISUALIZATION FOR DIFFERENTIAL EXPRESSION FINISHED")
    }

    setwd("..")

    #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #                                                                                       ## LASSO Regression ###
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    print("PROCESSING LASSO")

    #Lasso

    if (!isFALSE(lasso)) {
        if(lasso <= 0.0 || lasso > 1.0 ) {
            stop("\n- The input for argument lasso denotes hyperparameter alpha. This value must be set to 0.0 < x < 1.0 for Elastic Net (0.5 is default) or to 1.0 for LASSO regression. Re-run the pipeline again with correct lasso input or remove lasso all together.\n")
        }


        # Length of each group
        len <- as.numeric(table(group1))
        test.train <- unique(len >= 19)
        too.few <- unique(len < 9)


        # Stop Lasso if too few samples
        if (TRUE %in% too.few) {
            stop("\n- LASSO cannot be performed, too few samples per group, minimum is 10!\n")
        }


        dir.create("LASSOResults")
        setwd("LASSOResults/")

        group1.LASSO <- group1
        seeds <- sample(1:1000, 10)
        LASSO.res <- list()

        cat("Cross-validation for grouped multinomial LASSO is running with 10 random seeds, this will take some minutes...")

        if(FALSE %in% test.train) {
            if (databatch1 == TRUE) {
                if (length(levels(as.factor(group1.LASSO))) > 2) {
                    for (idx in 1:length(seeds)) {
                        LR <- LASSOFeature(seeds[[idx]], data1.batch, group1.LASSO, lasso, FALSE ,TRUE)
                        LASSO.res[[idx]] <-  LR
                    }
                } else {
                    for (idx in 1:length(seeds)) {
                        LR <- LASSOFeature(seeds[[idx]], data1.batch, group1.LASSO, lasso, FALSE, FALSE)
                        LASSO.res[[idx]] <-  LR
                    }
                }
            } else {
                if (length(levels(as.factor(group1.LASSO))) > 2) {
                    for (idx in 1:length(seeds)) {
                        LR <- LASSOFeature(seeds[[idx]], data1, group1.LASSO, lasso, FALSE ,TRUE)
                        LASSO.res[[idx]] <-  LR
                    }
                } else {
                    for (idx in 1:length(seeds)) {
                        LR <- LASSOFeature(seeds[[idx]], data1, group1.LASSO, lasso, FALSE ,FALSE)
                        LASSO.res[[idx]] <-  LR
                    }
                }
            }
        } else {
            if (databatch1 == TRUE) {
                if (length(levels(as.factor(group1.LASSO))) > 2) {
                    for (idx in 1:length(seeds)) {
                        LR <- LASSOFeature(seeds[[idx]], data1.batch, group1.LASSO, lasso, TRUE ,TRUE)
                        LASSO.res[[idx]] <-  LR
                    }
                } else {
                    for (idx in 1:length(seeds)) {
                        LR <- LASSOFeature(seeds[[idx]], data1.batch, group1.LASSO, lasso, TRUE, FALSE)
                        LASSO.res[[idx]] <-  LR
                    }
                }
            } else {
                if (length(levels(as.factor(group1.LASSO))) > 2) {
                    for (idx in 1:length(seeds)) {
                        LR <- LASSOFeature(seeds[[idx]], data1, group1.LASSO, lasso, TRUE ,TRUE)
                        LASSO.res[[idx]] <-  LR
                    }
                } else {
                    for (idx in 1:length(seeds)) {
                        LR <- LASSOFeature(seeds[[idx]], data1, group1.LASSO, lasso, TRUE ,FALSE)
                        LASSO.res[[idx]] <-  LR
                    }
                }
            }
        }


        # Extract results of 10 runs - Write out and plot results
        VarsSelect <- Reduce(intersect, lapply(LASSO.res, '[[', 1))

        if (length(VarsSelect) < 2) {
            stop("\n- There is no overlap in 10 elastic net runs. If you ran LASSO (lasso was et to 1.0) you can try and relax alpha and perform elastic net instead (0.0 < lasso < 1.0). Otherwise you data may have to high of a noise ratio to sample size, LASSO should not be performed.\n")
        }


        VarsSelect <- data.frame(VarsSelect[-1])
        colnames(VarsSelect) <- c("LASSO.Var.Select")


        # Write out LASSO/EN results
        write.table(VarsSelect, paste0(prefix,"_LASSO.txt"), sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


        # Consensus DEA and LASSO
        consensus <- DEA1.out[DEA1.out$name %in% VarsSelect$LASSO.Var.Select,]

        if (nrow(consensus) > 0) {
            write.table(consensus, paste0(prefix,"_DEA_LASSO_Consensus.txt"), sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
            pdf(paste0(prefix, "_overlap_DEAA_LASSO_EN.pdf"), height=8, width=12)
            if (length(levels(group1)) == 2) {
                venn <- venn.diagram(list(A=unique(as.character(DEA1.out[DEA1.out$dir =="up",]$name)), B=unique(as.character(DEA1.out[DEA1.out$dir =="down",]$name)), C=as.character(VarsSelect$LASSO.Var.Select)), category.names = c("D(EA) Analysis Up", "D(EA) Analysis Down", "LASSO/EN Regression"), filename=NULL, lwd = 0.7, cat.pos=0, sub.cex = 2, cat.cex= 1.5, cex=1.5, fill=viridis(3, begin = 0.2, end = 0.8, option="cividis"))

            } else {
                venn <- venn.diagram(list(A=unique(as.character(DEA1.out$name)), B=as.character(VarsSelect$LASSO.Var.Select)), category.names = c("D(EA) Analysis All", "LASSO/EN Regression"), filename=NULL, lwd = 0.7, cat.pos=0, sub.cex = 2, cat.cex= 1.5, cex=1.5, fill=viridis(2, begin = 0.2, end = 0.8, option="cividis"))

            }
            grid.draw(venn)
            dev.off()
        } else {
            cat("There is no consensus between LASSO regression and DEA/DAA.")
        }



        # Cross Validation errors
        LassoRun <- paste0(rep("Run", 10), 1:10)
        CrossValErrormean <- round(unlist(lapply(LASSO.res, '[[', 2)), digits = 4)
        cat(paste0("\nThe average leave-one-out cross validation error for LASSO/elastic-net was: ", mean(CrossValErrormean), "% and the higest error returned from any of the 10 runs was: ", max(CrossValErrormean),"%. Generally the cross validation error should be low ~ 5.0 %, as large errors may indicate a poor model and/or very heterogeneous data. On the other hand, an error of 0 might indicate over-fitting. See CAMPP manual for specifics.\n\n"))
        pCVEM <- data.frame(cbind(CrossValErrormean, LassoRun))
        pCVEM <- ggplot(data=pCVEM, aes(x=LassoRun, y=CrossValErrormean)) + geom_bar(aes(fill = as.factor(LassoRun)), stat="identity") + theme_minimal() + scale_x_discrete(limits=c(LassoRun)) + scale_fill_viridis(begin = 0.0, end = 0.0, discrete=TRUE, option="cividis" ) + theme(legend.position="none") + ylab("CrossValErrormean in %") + theme(axis.text = element_text(size=14), axis.title = element_text(size=16))
        ggsave(paste0(prefix, "_CrossValidationPlot.pdf"), plot = pCVEM)




        # Area under the curve AUC
        if(TRUE %in% test.train) {

            ll <- list()
            llev <- levels(as.factor(group1.LASSO))

            for (idx in 1:length(llev)) {
                pos <- which(group1.LASSO == as.character(llev[idx]))
                ll[[idx]] <- pos
            }

            my.samp <- unlist(lapply(ll, function(x) sample(x, ceiling((length(x)/4)))))


            if (databatch1 == TRUE) {
                LASSO.data1 <- data1.batch
            } else {
                LASSO.data1 <- data1
            }


            testD <- data.frame(t(LASSO.data1[rownames(LASSO.data1) %in% as.character(VarsSelect$LASSO.Var.Select), my.samp]))
            testG <- as.integer(group1.LASSO[my.samp])

            trainD <- data.frame(t(LASSO.data1[rownames(LASSO.data1) %in% as.character(VarsSelect$LASSO.Var.Select), -my.samp]))
            trainG <- as.integer(group1.LASSO[-my.samp])

            mn.net <- nnet::multinom(trainG ~ ., data=trainD)
            mn.pred <- predict(mn.net, newdata=testD, type="prob")
            roc.res <- multiclass.roc(testG, mn.pred)
            roc.res <- data.frame(round(as.numeric(sub(".*: ", "", roc.res$auc)), digits = 2))
            colnames(roc.res) <- "AUC"
            cat(paste0("Are under the curve (AUC) for variables selected from 10 LASSO/EN runs was: ", roc.res$AUC))
            write.table(roc.res, paste0(prefix,"_AUC.txt"), row.names=FALSE, col.names = TRUE, quote = FALSE)
        }


        setwd("..")
        try(rm(VarsSelect, LR, venn, CrossValErrormean, mn.net, mn.pred, roc.res), silent=T)

    }

    print("LASSO PART FINISHED")

    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ### Plotting Results Heatmap ###
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    if (!isFALSE(plot.heatmap)) {

        print("PRINTING HEATMAP FOR FIRST DATASET")

        if (databatch1 == TRUE){
            hm <- data1.batch
            name <- "_batchcorr"
        } else {
            hm <- data1
            name <- ""
        }

        if (plot.heatmap == "ALL") {
            stop("\nOption ALL is not allowed for heatmap, too heavy! Pick either 'DEA', 'LASSO', 'EN' or 'Consensus'.\n")
        } else if (plot.heatmap %in% c("DEA")) {
            signif_samples <- rbind(head(DE.out[order(DE.out$logFC),],heatmap.size/2),tail(DE.out[order(DE.out$logFC),],heatmap.size/2))
            signif_samples <- signif_samples$name
            hm <- data1[rownames(data1) %in% signif_samples,]
            cat("\n- DEA was selected. Normalized and voom transformed data will be used for visualization.\n")
        } else {
            if(!is.null(lasso)) {
                if (plot.heatmap == "Consensus") {
                    hm <- hm[rownames(hm) %in% as.character(consensus$name),]
                }
                if (plot.heatmap %in% c("EN", "LASSO")) {
                    hm <- hm[rownames(hm) %in% as.character(VarsSelect$LASSO.Var.Select),]
                }
            } else {
                stop("You have specified argument plot.heatmap which will produce a heatmap with results from either DA/DE analysis, LASSO/Elastic-Net Regression or the Consensus of these. Input to argument plot.heatmap must be a string specifying which results to plot, options are: DA, DE, LASSO, EN, Consensus")
            }
        }

        # Heatmap colors in blue
        hm.gradient <- viridis(300, option="cividis")
        range <- c(round(min(hm)), round(max(hm)))

        # Add HUGO IDs
        name <- rownames(hm)
        rownames(hm) <- NULL
        hm <- cbind(name,hm)
        hm <- AddGeneName(hm,ensembl.version)
        hm <- subset(hm, HUGO_ID != "")
        rownames(hm) <- hm$HUGO_ID
        hm <- subset(hm, select=-c(Ensembl_ID,HUGO_ID))

        # Heatmap as pdf
        MakeHeatmap(hm, hm.gradient, group1, prefix, range)
        rm(hm)

        print("HEATMAP WAS SUCCESFULLY PRINTED")

        if (!is.null(data2)){
            print("PRINTING HEATMAP FOR SECOND DATASET")

            if (databatch2 == TRUE){
                hm2 <- data2.batch
                name <- "_batchcorr"
            } else {
                hm2 <- data2
                name <- ""
            }

            if (plot.heatmap %in% c("DEA")) {
                signif_samples2 <- rbind(head(sDE.out[order(sDE.out$logFC),],heatmap.size/2),tail(sDE.out[order(sDE.out$logFC),],heatmap.size/2))
                signif_samples2 <- signif_samples2$name
                hm2 <- data2[rownames(data2) %in% signif_samples2,]
                cat("\n- DEA was selected. Normalized and voom transformed data will be used for visualization.\n")
            } else {
                stop("A heatmap for the second dataset can only be generated for DEA")
            }

            # Add HUGO IDs
            name <- rownames(hm2)
            rownames(hm2) <- NULL
            hm2 <- cbind(name,hm2)
            hm2 <- AddGeneName(hm2,ensembl.version)
            hm2 <- subset(hm2, HUGO_ID != "")
            rownames(hm2) <- hm2$HUGO_ID
            hm2 <- subset(hm2, select=-c(Ensembl_ID,HUGO_ID))

            range <- c(round(min(hm2)), round(max(hm2)))

            # Heatmap as pdf
            MakeHeatmap(hm2, hm.gradient, group2, paste0(prefix,"_second"), range)
            rm(hm2)

            print("HEATMAP WAS SUCCESFULLY PRINTED")
        }

    } else {
        cat("\n- No heatmap requested.\n")
    }


    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ### CORRELATION ANALYSIS ###
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    print("PROCESSING CORRELATION ANALYSIS")

    if (!is.null(data2) & !isFALSE(corrby)) {

        dir.create("CorrelationResults")
        setwd("CorrelationResults/")


        retainedSamp <- intersect(colnames(data1), colnames(data2))

        if (corrby == "ALL") {
            retainedvar <- intersect(rownames(data1), rownames(data2))
        } else if (corrby == "DEA") {
            retainedvar <- intersect(res.DE.names, rownames(data2))
        } else {
            if(!is.null(lasso)) {
                if (survival == "Consensus" ) {
                    retainedvar <- intersect(as.character(consensus$name), rownames(data2))
                } else {
                    retainedvar <- intersect(as.character(VarsSelect$LASSO.Var.Select), rownames(data2))
                }
            } else {
                stop("You have chosen to perform correlation analysis with results from LASSO/EN BUT this analysis has not been performed. Please re-run pipeline with parameter lasso (see Manual or help (-h))")
            }
        }


        if (databatch1 == TRUE) {
            data1.corr <- data1.batch[rownames(data1.batch) %in% retainedvar, colnames(data1.batch) %in% retainedSamp]
        } else {
            data1.corr <- data1[rownames(data1) %in% retainedvar, colnames(data1) %in% retainedSamp]
        }

        if (databatch2 == TRUE) {
            corr.out <- data2.batch[rownames(data2.batch) %in% retainedvar, colnames(data2.batch) %in% retainedSamp]
        } else {
            corr.out <- data2[rownames(data2) %in% retainedvar, colnames(data2) %in% retainedSamp]
        }

        # Perform correlation analysis and generate overall correlation plot
        res.corr <- CorrAnalysis(data1.corr, corr.out, prefix)


        # print out significant hits in Excel
        res.corr$sig.corr <- ifelse(res.corr$fdr <= FDR1, "yes", "no")

        write.table(res.corr, paste0(prefix,"_corr.txt"), sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


        # Individual correlation plots
        corr.features <- as.character(res.corr[which(res.corr$sig.corr == "yes"),]$name)
        if (length(corr.features) > 1) {
            CorrelationPlots(data1.corr, corr.out, corr.features, prefix)
        } else {
            cat("\n- No significant correlations, no plotting.\n")
        }
        setwd("..")
        try(rm(res.corr,corr.out,data1.corr), silent=T)
    }


    print("CORRELATION ANALYSIS PART FINISHED")


    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    ### SURVIVAL ANALYSIS ###
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    print("PROCESSING SURVIVAL ANALYSIS")



    ### SETTING UP DATA ###


    if (!isFALSE(survival)) {

        if (length(intersect(colnames(metadata1), survival.metadata)) < 3) {
            stop("Survival analysis requested, but one or more columns; 'survival', 'outcome', 'outcome.time' are missing from metadata file!")
        }


        dir.create("SurvivalResults")
        setwd("SurvivalResults/")

        # Setting up dataset for survival analysis, matching samples either batch or no batch. Only DA/DE features are used.
        metadata1.surv <- metadata1[which(metadata1$survival == 1),]
        samples <- metadata1.surv$ids

        if (databatch1 == TRUE) {
            data1.surv <- data1.batch[,colnames(data1.batch) %in% samples]
        } else {
            data1.surv <- data1[,colnames(data1) %in% samples]
        }



        if (survival == "ALL") {
            cat("\nYou have chosen to use all varibles for survival analysis, this is NOT advisable! See options for survival in with the help argument.\n")
        } else if (survival == "DEA") {
            data1.surv <- data1.surv[rownames(data1.surv) %in% res.DE.names,]
        } else {
            if(!is.null(lasso)) {
                if (survival == "Consensus" ) {
                    data1.surv <- data1.surv[rownames(data1.surv) %in% as.character(consensus$name),]
                } else {
                    data1.surv <- data1.surv[rownames(data1.surv) %in% as.character(VarsSelect$LASSO.Var.Select),]
                }
            } else {
                stop("You have chosen to perform survival analysis with results from LASSO/EN but this analysis has not been performed. Please re-run pipeline with parameter lasso (see Manual)")
            }
        }


        # If user has specified covariates for survival analysis, extract these.
        if(is.null(covarS)) {
            surv_object <- data.frame(t(data1.surv), as.numeric(metadata1.surv$outcome.time), metadata1.surv$outcome)
            colnames(surv_object) <- c(rownames(data1.surv), "outcome.time", "outcome")
        } else {
            my.covars <- metadata1.surv[,colnames(metadata1.surv) %in% covarS]
            surv_object <- data.frame(t(data1.surv), as.numeric(metadata1.surv$outcome.time), metadata1.surv$outcome, my.covars)
            colnames(surv_object) <- c(rownames(data1.surv), "outcome.time", "outcome", covarS)
        }


        # datadist format for cph function with cubic splines.
        features <- rownames(data1.surv)
        dd <- datadist(surv_object)
        options(datadist=dd)

        mylist <- list(features, dd, surv_object)


        ### User Specified Covariates ###


        # Evaluating class of covariate, continious or categorical.
        covarS.original <- covarS

        if(!is.null(covarS)) {
            for (i in 1:length(covarS)) {
                if(class(eval(parse(text=paste0("surv_object$", covarS[i])))) %in% c("integer", "numeric")) {
                    covarS[i] <- paste0("rcs(", covarS[i], ")")
                }
            }

            covarS <- paste(covarS,collapse="+")
        }


        ### Liniarity of Continious Covariates ###


        # Survival models with or without user specified covariates.
        covariate_linearity <- list()

        if(!is.null(covarS)) {
            for (f in features) {
                acall <- parse(text = paste0("result <- cph(Surv(outcome.time, outcome) ~ rcs(", as.character(f),") + ", covarS, ", data = surv_object, x=TRUE,y=TRUE)"))
                eval(acall)
                covariate_linearity[[as.character(f)]] <- result
            }
        } else {
            for (f in features) {
                acall <- parse(text = paste0("result <- cph(Surv(outcome.time, outcome) ~ rcs(", as.character(f),"), data = surv_object, x=TRUE,y=TRUE)"))
                eval(acall)
                covariate_linearity[[as.character(f)]] <- result
            }
        }


        # Check for liniarity violation:
        check1 <- which(as.numeric(unlist(lapply(covariate_linearity, function(x) length(x)))) == 1)

        if (length(check1) > 0) {
            surv_object <- surv_object[,-check1]
            covariate_linearity <- covariate_linearity[-check1]
        }

        # Check for error, other.
        check2 <- lapply(covariate_linearity, function(i) tryCatch({ anova(i) }, error=identity))
        check2 <- which(as.character(vapply(check2, is, logical(1), "error")) == TRUE)


        if (length(check2) > 0) {
            surv_object <- surv_object[,-check2]
            covariate_linearity <- covariate_linearity[-check2]
        }


        # Remove checks to save memory
        rm(check1, check2)

        # Re-assign features
        features <- colnames(surv_object)

        # Anova comparing linear covariate and non-linear covariate.
        covariate_linearity  <- lapply(covariate_linearity, function(x) anova(x))

        # Unlisting and extracting pvalues for non-linear fit.
        covariate_linearity <- data.frame(do.call(rbind, lapply(covariate_linearity, function(x) x[grep("Nonlinear|NONLINEAR", rownames(x)),3])))

        # FDR correction for multiple testing
        covariate_linearity <- apply(covariate_linearity, 2, function(x) p.adjust(x, method = "fdr", n=nrow(covariate_linearity)))

        # Colnames with and without user-specified covariates.
        if(is.null(covarS)) {
            colnames(covariate_linearity) <- c("Feature")
        } else {
            covar.names <- unlist(strsplit(covarS, split="[+]"))
            if(length(grep("rcs", covar.names)) > 0) {
                covar.names <- c("Feature", covar.names[grep("rcs", covar.names)], "Total")
                covar.names <- gsub("rcs\\(|\\)", "", covar.names)
                colnames(covariate_linearity) <- covar.names
            } else {
                colnames(covariate_linearity) <- c("Feature")
            }
        }


        # Extracting p-values and filtering for significance.
        covariate_nonlinear <- colnames(covariate_linearity)[apply(covariate_linearity, 2, function(x) any(x < 0.05))]

        # Updating covariates with cubic splines.
        if (length(covariate_nonlinear) > 0) {
            cat(paste0("\n- The following continious covariate(s) may be violating the assumption of linearity:  ", covariate_nonlinear,".\nCubic splines will be added.\n"))
            if(!is.null(covarS.original)) {
                for (i in 1:length(covarS.original)) {
                    if (covarS.original[i] %in% covariate_nonlinear) {
                        covarS.original[i] <- paste0("rcs(", covarS.original[i], ")")
                    }
                }
            }
        }

        rm(covariate_linearity)


        ### Testing Cox Proportional Hazard Assumption ###


        # Picking model based on result of linearity check.
        pha_check <- list()


        for (f in features) {
            if ("feature" %in% covariate_nonlinear) {
                if(is.null(covarS.original)) {
                    acall <- parse(text = paste0("result <- data.frame(cox.zph(cph(Surv(outcome.time, outcome) ~ rcs(", as.character(f), "), data = surv_object, x=TRUE,y=TRUE))$table[,3])"))
                    eval(acall)
                    pha_check[[as.character(f)]] <- result
                } else {
                    acall <- parse(text = paste0("result <- data.frame(cox.zph(cph(Surv(outcome.time, outcome) ~ rcs(", as.character(f),") + ", covarS.original, ", data = surv_object, x=TRUE,y=TRUE))$table[,3])"))
                    eval(acall)
                    pha_check[[as.character(f)]] <- result
                }
            } else {
                if(is.null(covarS.original)) {
                    acall <- parse(text = paste0("result <- data.frame(cox.zph(cph(Surv(outcome.time, outcome) ~ ", as.character(f), ", data = surv_object, x=TRUE,y=TRUE))$table[,3])"))
                    eval(acall)
                    pha_check[[as.character(f)]] <- result

                } else {
                    acall <- parse(text = paste0("result <- data.frame(cox.zph(cph(Surv(outcome.time, outcome) ~ ", as.character(f), "+ ", covarS.original, ", data = surv_object, x=TRUE,y=TRUE))$table[,3])"))
                    eval(acall)
                    pha_check[[as.character(f)]] <- result
                }
            }
        }




        # Filtering p-values for significance.
        if(length(pha_check[1]) == 1) {
            pha.fail.test <- unlist(lapply(pha_check, function(x) any(x < 0.05)))
            pha.fail.test <- names(pha.fail.test)[pha.fail.test == TRUE]
        } else {
            pha.fail.test <- as.character(unlist(lapply(pha_check, function(x) rownames(x)[apply(x, 1, function(u) any(u < 0.05))])))
            if ("GLOBAL" %in% pha.fail.test) {
                pha.fail.test <- pha.fail.test[-which(pha.fail.test == "GLOBAL")]
            }
        }



        cat("\nWARNING: The following features and/or covariates failed the test of proportional hazard: ", pha.fail.test, "\nIF the covariates that failed are categorical you may use strata by re-running the pipline adding argument stratify followed by the names of the categorical covariates to stratify (if multiple then separate by comma). \nN.B, this pipeline does not handle continuous variables that violate the proportional hazard assumption, if any of these failed PH test, the hazard ratios of these should NOT be evaluated.\n")



        # User stratification of categorical covariates which violate the proportional hazard assumption.
        if(!is.null(stratify)) {
            for (i in 1:length(covarS.original)) {
                if (covarS.original[i] %in% stratify) {
                    covarS.original[i] <- paste0("strat(", covarS.original[i], ")")
                }
            }
        }


        ### Survival Analysis with updated covariates ###
        # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        # Survival analysis, cox-regression


        if(!is.null(covarS)) {
            features.surv <- features[-c((length(features)-length(covariates)):length(features))]
        } else {
            features.surv <- features[-c((length(features)-1):length(features))]
        }

        # Picking model based on result of linearity check and proportional hazard text.
        survival.results <- list()



        for (f in features.surv) {
            if ("feature" %in% covariate_nonlinear) {
                if(is.null(covarS.original)) {
                    acall <- parse(text = paste0("result <- cph(Surv(outcome.time, outcome) ~ rcs(", as.character(f),"), data = surv_object, x=TRUE,y=TRUE)"))
                    eval(acall)
                    survival.results[[as.character(f)]] <- result
                } else {
                    acall <- parse(text = paste0("result <- cph(Surv(outcome.time, outcome) ~ rcs(", as.character(f),") +", covarS.original, ", data = surv_object, x=TRUE,y=TRUE)"))
                    eval(acall)
                    survival.results[[as.character(f)]] <- result
                }
            } else {
                if(is.null(covarS.original)) {
                    acall <- parse(text = paste0("result <- cph(Surv(outcome.time, outcome) ~", as.character(f), ", data = surv_object, x=TRUE,y=TRUE)"))
                    eval(acall)
                    survival.results[[as.character(f)]] <- result

                } else {
                    acall <- parse(text = paste0("result <- cph(Surv(outcome.time, outcome) ~", as.character(f), " +", covarS.original, ", data = surv_object, x=TRUE,y=TRUE)"))
                    eval(acall)
                    survival.results[[as.character(f)]] <- result
                }
            }
        }


        # Setting up data and writing out excel file with results and making HR stemplot
        survival.data1 <- SurvivalCOX(survival.results, prefix, surv.plot)

        write.table(survival.data1, paste0(prefix,"_survival.txt"), sep = "\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

        setwd("..")
        rm(dd, pha_check, pha.fail.test, survival.results, features.surv, survival.data1)
    }

    print("SURVIVAL ANALYSIS PART FINISHED")


    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #                                                               ## Weighed Gene Co-Expression Network Analysis ###
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    print("PROCESSING WGCNA")

    if (!isFALSE(WGCNA)) {

        dir.create("WGCNAResults")
        setwd("WGCNAResults/")

        if (databatch1 == TRUE){
            data1.WGCNA <- t(data1.batch)
        } else {
            data1.WGCNA <- t(data1)
        }

        if (WGCNA == "DEA") {
            data1.WGCNA <- data1.WGCNA[,colnames(data1.WGCNA) %in% res.DEA1.names]
        }

        dir.create("WGCNAPlots")
        setwd("WGCNAPlots/")
        # Check data
        gsg <- goodSamplesGenes(data1.WGCNA)
        cat(paste0("- Data set is OK for WGCNA - ", gsg$allOK,".\n"))
        WGCNAres <- WGCNAAnalysis(data1.WGCNA, cutoff.WGCNA, prefix)
        setwd("..")

        if (WGCNA %in% c("DA", "DE")) {
            logFCs <- DEA1.out
            logFCs$gene <- DEA1.out$name
            WGCNAres <- lapply(WGCNAres, function(x) merge(x, logFCs, by = "gene"))
            WGCNAres <- lapply(WGCNAres, function(x) x[with(x, order(comparison, Module)),])
        }

        mod.names <- names(WGCNAres)
        for (idx in 1:length(WGCNAres)) {
            write.table(WGCNAres[[idx]], paste0(mod.names[idx],"_moduleRes.txt"), sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
        }

        setwd("..")
        rm(WGCNAres)

    }

    gc()
    print("WGCNA PART FINISHED")


    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    #                                                               ## Protein-Protein and MiRNA-Gene Network Analysis ###
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    print("PROCESSING PROTEIN-PROTEIN and MiRNA-GENE NETWORK ANALYSIS")

    if (!is.null(PPI)) {

        DB <- DownloadPPInt(PPI[[1]])

        dir.create("InteractionResults")
        setwd("InteractionResults/")
        dir.create("InteractionPlots")

        DEA1.out$name <- gsub("_", "-", DEA1.out$name)
        PPIs <- GetDeaPPInt(DB, DEA1.out)

        if (!is.null(GmiRI)) {

            DEA2.out$name <- gsub("_", "-", DEA2.out$name)

            cat("\nSearching for miRNA-gene interaction - this may take up to 10 min!\n")
            GmiRIs <- GetGeneMiRNAInt(GmiRI[[2]], DEA2.out, DEA1.out)

            PPGmiRIs <- MatchPPGmiRInt(GmiRIs, PPIs)
            PPGmiRTrim <- TrimWriteInt(PPGmiRIs)

            setwd("InteractionPlots/")
            cat("\nPlotting set-wise networks - this may take up to 10 min!\n")
            PlotInt(PPGmiRTrim)

            setwd("../..")

        } else {

            PPTrim <- TrimWriteInt(PPIs)

            setwd("InteractionPlots/")
            cat("\nPlotting set-wise networks - this may take up to 10 min!\n")
            PlotInt(PPTrim)

            setwd("../..")
        }

    } else if (!is.null(GmiRI)) {

        dir.create("InteractionResults")
        setwd("InteractionResults/")
        dir.create("InteractionPlots")

        DEA1.out$name <- gsub("_", "-", DEA1.out$name)

        cat("\nSearching for miRNA-gene interaction - this may take up to 10 min!\n")
        GmiRIs <- GetGeneMiRNAInt(GmiRI[[2]], DEA1.out)
        GmiRTrim <- TrimWriteInt(GmiRIs)

        setwd("InteractionPlots/")
        cat("\nPlotting set-wise networks - this may take up to 10 min!\n")
        PlotInt(GmiRTrim)

        setwd("../..")

    } else {
        cat("\nNo Interaction Networks requested\n.")

        print("PROTEIN-PROTEIN and MiRNA-GENE NETWORK ANALYSIS PART FINISHED")
    }

    setwd("../")
    cat("\nCAMPP RUN DONE!\n")
}
