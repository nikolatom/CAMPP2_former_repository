#' @title CAMPP2 pipeline
#' @description CAMPP2 is a tool for quantitative data analysis
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
#' @import zeallot
#' @export
#' @seealso
#' @return CAMPP2 results
#' @examples \dontrun{
#' ...
#' }


runCampp2 <- function (data, metadata, sdata=NULL, smetadata=NULL, technology, groups, batches=NULL, datacheck=TRUE, standardize=NULL, transform=NULL, plotmds=NULL, plotheatmap=NULL, kmeans=NULL, signif=NULL, colors=NULL, prefix=NULL, corr=NULL, lasso=NULL, WGCNA=NULL, cutoffWGCNA=NULL, survival=NULL, covar=NULL, stratify=NULL, survplot=NULL, PPint=NULL, GenemiRint=NULL){

  ###parse input arguments and assign updated vales
  c(data,sdata,metadata,smetadata,technology,groups,group,sgroup,acall,ids,batches,databatch,sdatabatch,standardize,transform,datacheck,plotmds,kmeans,labels.kmeans,signif,logFC,FDR,slogFC,sFDR,colors,prefix,plotheatmap,corrby,lasso,WGCNA,cutoffWGCNA,survival,covarD,scovarD,covarS,stratify,survplot,PPI,GmiRI) %<-% parseArguments(data=data, metadata=metadata, sdata=sdata, smetadata=smetadata, groups=groups, technology=technology, prefix=prefix, batches=batches, datacheck=datacheck, standardize=standardize, transform=transform, plotmds=plotmds, plotheatmap=plotheatmap, kmeans=kmeans, signif=signif, colors=colors, corr=corr, lasso=lasso, WGCNA=WGCNA, cutoffWGCNA=cutoffWGCNA, survival=survival, covar=covar, stratify=stratify, survplot=survplot, PPint=PPint, GenemiRint=GenemiRint)

  ###this is an alternative to the previous step
  # arguments<-parseArguments2(data=data, metadata=metadata, sdata=sdata, smetadata=smetadata, groups=groups, technology=technology, prefix=prefix, batches=batches, datacheck=datacheck, standardize=standardize, transform=transform, plotmds=plotmds, plotheatmap=plotheatmap, kmeans=kmeans, signif=signif, colors=colors, corr=corr, lasso=lasso, WGCNA=WGCNA, cutoffWGCNA=cutoffWGCNA, survival=survival, covar=covar, stratify=stratify, survplot=survplot, PPint=PPint, GenemiRint=GenemiRint)

    ###should the variables be assigned in every function?
  # data<-arguments$data
  # technology<-arguments$technology
  # groups<-arguments$groups
  # group<-arguments$group
  # sgroup<-arguments$sgroup
  # acall<-arguments$acall
  # ids<-arguments$ids
  # batches<-arguments$batches
  # databatch<-arguments$databatch
  # sdatabatch<-arguments$sdatabatch
  # standardize<-arguments$standardize
  # transform<-arguments$transform
  # datacheck<-arguments$datacheck
  # plotmds<-arguments$plotmds
  # kmeans<-arguments$kmeans
  # labels.kmeans<-arguments$labels.kmeans
  # signif<-arguments$signif
  # logFC<-arguments$logFC
  # FDR<-arguments$FDR
  # slogFC<-arguments$slogFC
  # sFDR<-arguments$sFDR
  # colors<-arguments$colors
  # prefix<-arguments$prefix
  # plotheatmap<-arguments$plotheatmap
  # corrby<-arguments$corrby
  # lasso<-arguments$lasso
  # WGCNA<-arguments$WGCNA
  # cutoffWGCNA<-arguments$cutoffWGCNA
  # survival<-arguments$survival
  # covarD<-arguments$covarD
  # scovarD<-arguments$scovarD
  # covarS<-arguments$covarS
  # stratify<-arguments$stratify
  # survplot<-arguments$survplot
  # PPI<-arguments$PPI
  # GmiRI<-arguments$GmiRI


  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Create directory Results
  dir.create(prefix)
  setwd(paste0(prefix, "/"))


  ######CREATE A NEW BRANCH; ALSO ONE FOR PARSING THE ARGUMENTS!!!

  print("RUN TRANSFORMATION")
  # Check if data contains zeros and negative values.
  hasZeroD <- unique(as.vector(data == 0))
  hasNegD <- unique(as.vector(data < 0))

  if(transform[1] %in% c("log2", "log10", "logit")) {
    if (TRUE %in% hasNegD) {
      stop("\n- Data contains negative values and cannot be log transformed. Re-run command WITHOUT argument transform  or alternatively if using two datasets, specify 'none' as the transforminput for the dataset with negative values, e.g. 'none,log2' or 'log2,none'.\n")
    } else {
      if (TRUE %in% hasZeroD) {
        data.original <- data
        data <- ReplaceZero(data, group)
      }
    }
  }


  if (!is.null(sdata)){
    hasZeroS <- unique(as.vector(sdata == 0))
    hasNegS <- unique(as.vector(sdata < 0))
  }

  if(!is.null(sdata) & transform[2] %in% c("log2", "log10", "logit")) {
    if (TRUE %in% hasNegS) {
      stop("\n- Second dataset contains negative values and cannot be log transformed. Re-run command WITHOUT argument transform  or alternatively if using two datasets, specify 'none' as the transforminput for the dataset with negative values, e.g. 'none,log2' or 'log2,none'.\n")
    } else {
      if (TRUE %in% hasZeroS) {
        sdata.original <- sdata
        sdata <- ReplaceZero(sdata, group)
      }
    }
  }


  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print("TEST GROUP")
print(group)

  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  #                                                                         ## Normalization, Filtering and Transformation ###
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  NB <- " N.B This pipeline does not handle background correction of single-channel intensity data or within array normalization two-color intensity data. See limma manual section on array normalization for more on this. Data may be fully normalized with limma (R/Rstudio) or another software and the pipeline re-run."


  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Dataset

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



  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ### BATCH CORRECTION ###
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ###create a function "runDatabatchCorr" taking arguments (databatch, technology, sdatabatch)


  if (databatch == TRUE){
    if (length(batch) > 0) {
      design <-  model.matrix(~group)

      if (technology[1] == "seq") {
        data.batch <- ComBat(as.matrix(data$E), batch, design, par.prior=TRUE,prior.plots=FALSE)
      } else {
        data.batch <- ComBat(as.matrix(data), batch, design, par.prior=TRUE,prior.plots=FALSE)
      }

    } else {
      data.batch <- data
      cat("\n- No column names match specified batches for dataset.\n")
    }
  } else {
    cat("\n- No batch correction requested.\n")
  }




  if (sdatabatch == TRUE){
    if (length(sbatch) > 0) {
      sdesign <- model.matrix(~sgroup)

      if (technology[2] == "seq") {
        sdata.batch <- ComBat(as.matrix(sdata$E), sbatch, sdesign, par.prior=TRUE,prior.plots=FALSE)
      } else {
        sdata.batch <- ComBat(as.matrix(sdata), sbatch, sdesign, par.prior=TRUE,prior.plots=FALSE)
      }
    } else {
      sdata.batch <- sdata
      cat("\n- No column names match specified batches for second dataset. Continuing without batch correction.\n")
    }
  } else {
    cat("\n- No batch correction requested for second dataset.\n")
  }



  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ### Distributional Checks ###
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  if (datacheck == TRUE) {
    if (databatch == TRUE) {
      subset.data <- data.batch[sample(nrow(data.batch), 10),]
    } else {
      if (technology[1] == "seq") {
        subset.data <- data$E[sample(nrow(data$E), 10),]
      } else {
        subset.data <- data[sample(nrow(data), 10),]
      }
    }

    list.of.dists <- FitDistributions(subset.data)

    dir.create("DataChecks")
    setwd("DataChecks/")
    PlotDistributions(subset.data, list.of.dists)
    setwd("..")

    rm(subset.data, list.of.dists)
  }



  if (datacheck == TRUE & !is.null(sdata)) {
    if (sdatabatch == TRUE) {
      subset.data <- sdata.batch[sample(nrow(sdata.batch), 10),]
    } else {

      if (technology[2] == "seq") {
        subset.data <- sdata$E[sample(nrow(sdata$E), 10),]
      } else {
        subset.data <- sdata[sample(nrow(sdata), 10),]
      }
    }

    list.of.dists <- FitDistributions(subset.data)

    dir.create("SecondDataChecks")
    setwd("SecondDataChecks/")
    PlotDistributions(subset.data, list.of.dists)
    setwd("..")

    rm(subset.data, list.of.dists)
  }



  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------






  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ### PRELIMINARY MDS PLOT ###
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



  MDScolors <- gsub(pattern = "FF", replacement = "", x = colors)


  if (plotmds == TRUE && databatch == TRUE){
    mdsplot <- MDSPlot(data.batch, group, ids, MDScolors)
    ggsave(paste0(prefix, "_MDSplot_batchcorr.pdf"), plot = mdsplot, dpi = 300, width = 8, height = 8)

  } else if (plotmds == TRUE && databatch == FALSE){
    if (technology[1] == "seq") {
      mdsplot <- MDSPlot(data.frame(data$E), group, ids, MDScolors)
    } else {
      mdsplot <- MDSPlot(data, group, ids, MDScolors)
    }
    ggsave(paste0(prefix, "_MDSplot.pdf"), plot = mdsplot, dpi = 300, width = 8, height = 8)

    rm(mdsplot)

  } else {
    cat("\n- No preliminary plot requested.\n")
  }



  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ### Kmeans ###
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





  if (kmeans == TRUE) {

    if (technology[1] == "seq") {
      k.data <- data$E
    } else {
      k.data <- data
    }

    # Number of sample sets to generate
    nsets <- 1:ceiling(nrow(k.data)/1000)
    if(length(nsets) > 10) {
      nsets <- 1:10
    }

    # Number of variables in each sample set
    setsize <- nrow(k.data)
    if (setsize > 2000) {
      setsize <- 2000
    }


    # Number of kmeans to try
    if(ncol(k.data) <= 100) {
      nks <- 2:6
    } else if (ncol(k.data) > 100 && ncol(k.data) <= 500) {
      nks <- 2:11
    } else {
      nks <- 2:16
    }


    cat(paste0("Based on size of dataset, ", length(nsets), " sample sets will be generated of size ", setsize, " and ", length(nks), " clusters will be tested - N.B This may take up to 10 min!\nRunning......"))

    dir.create("KmeansResults")
    setwd("KmeansResults/")


    list.of.dfs <- list()

    if (databatch == TRUE) {
      for (idx in 1:length(nsets)) {
        df <- t(data.batch[sample(nrow(data.batch), setsize), ])
        list.of.dfs[[idx]] <- df
      }
      Kmeans.list <- lapply(list.of.dfs, function(x) EstimateKmeans(x, nsets))
      Kmeans.Out <- PlotKmeans(data.batch, Kmeans.list, nks, labels.kmeans, prefix)
    } else {
      for (idx in 1:length(nsets)) {
        df <- t(k.data[sample(nrow(k.data), setsize), ])
        list.of.dfs[[idx]] <- df
      }
      Kmeans.list <- lapply(list.of.dfs, function(x) EstimateKmeans(x, nsets))
      Kmeans.Out <- PlotKmeans(k.data, Kmeans.list, nks, labels.kmeans, prefix)
    }


    out <- cbind(metadata, Kmeans.Out)

    write.table(out, paste0(prefix,"_Metadata_Kmeans.txt"), sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

    setwd("..")

    rm(k.data, Kmeans.list, Kmeans.Out, out)
  }






  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # CLEAN UP AND SOURCE FUNCTIONS FROM THE FUNCTIONS SCRIPT - PART 2
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  #rm(SplitList, ReadMyFile, ReplaceNAs, ReplaceZero, NormalizeData, FitDistributions, PlotDistributions, MDSPlot, EstimateKmeans)
  gc(full = TRUE)



  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ### DIFFERENTIAL EXPRESSION ANALYSIS ###
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



  # Differential Expression Analysis with Limma


  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # First dataset

  dir.create("DEAAResults")
  setwd("DEAAResults/")


  # Make design matrix
  if (databatch == "FALSE") {
    design.str <- "model.matrix(~0+group"
    out.name <- "_DE"
  } else if (length(batch) > 0) {
    design.str <- "model.matrix(~0+group+batch"
    out.name <- "_databatch_DE"
  } else {
    stop("Batch correction selected but no batches column found!")
  }

  if(is.null(covarD)) {
    design <- eval(parse(text=paste0(design.str, ")")))
  } else {
    if (length(covarD) == 1) {
      df <- data.frame(metadata[,colnames(metadata) %in% covarD])
      colnames(df) <- covarD
    } else {
      df <- metadata[,colnames(metadata) %in% covarD]
    }

    s <- lapply(split(as.matrix(df), col(df)), factor)
    my.names <- paste0("", colnames(df))
    list2env(setNames(s, my.names), envir=.GlobalEnv)
    my.names <- paste0(my.names, collapse = "+")
    design <- eval(parse(text=paste0(design.str,"+",my.names,")")))
    rm(s)
  }



  # Making group contrasts
  combinations <- data.frame(t(combn(paste0("group", levels(group)), 2)))
  combinations$contr <- apply(combinations[,colnames(combinations)], 1, paste, collapse = "-")
  contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list(as.character(combinations$contr)),levels=list(design))))


  # Apply DE_limma function to all comparisons
  res.DE <- DAFeatureApply(contrast.matrix, data, design, logFC, FDR, NULL, FALSE)


  # Write results out as excel file
  if (!is.null(res.DE)) {
    #DE.out <- ExcelOutput(res.DE, paste0(prefix, out.name))
    DE.out <- TextOutput(res.DE, paste0(prefix, out.name))
    rownames(DE.out) <- NULL
    res.DE.names <- unique(DE.out$name)
    rm(res.DE)
  } else {
    cat("No signficant DE/DA hits found. Check output file from differential expression analysis. Check your cut-off for differential expression analysis, it may be that these are too stringent.")
  }

  setwd("..")




  if (technology[1] == "seq") {
    cnames <- colnames(data$E)
    data <- data.frame(data$E)
    colnames(data) <- cnames
  }





  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Second dataset

  ##line 845 was originally "if(!is.null(arg.PPI) & !is.null(arg.GmiRI)) {"

  if(!is.null(sdata) & !is.null(smetadata)) {
    setwd("DEAAResults/")

    combinations <- data.frame(t(combn(paste0("sgroup", levels(sgroup)), 2)))
    combinations$contr <- apply(combinations[,colnames(combinations)], 1, paste, collapse = "-")


    if (sdatabatch == FALSE) {
      sdesign.str <- "model.matrix(~0+sgroup"
      out.name <- "_second_DE"
    } else if (length(sbatch) > 0) {
      sdesign.str <- "model.matrix(~0+sgroup+sbatch"
      out.name <- "_second_databatch_DE"
    } else {
      stop("Batch correction selected but no batches column found!")
    }

    if(is.null(scovarD)) {
      sdesign <- eval(parse(text=paste0(sdesign.str, ")")))
    } else {
      if (length(scovarD) == 1) {
        df <- data.frame(smetadata[,colnames(smetadata) %in% scovarD])
        colnames(df) <- scovarD
      } else {
        df <- smetadata[,colnames(smetadata) %in% scovarD]
      }

      s <- lapply(split(as.matrix(df), col(df)), factor)
      my.names <- paste0("s", colnames(df))
      list2env(setNames(s, my.names), envir=.GlobalEnv)
      my.names <- paste0(my.names, collapse = "+")
      sdesign <- eval(parse(text=paste0(sdesign.str,"+",my.names,")")))
      rm(s)
    }

    # Making group contrasts
    contrast.matrix <- eval(as.call(c(as.symbol("makeContrasts"),as.list(as.character(combinations$contr)),levels=list(sdesign))))

    # Apply DE_limma function to all comparisons
    res.sDE <- DAFeatureApply(contrast.matrix, sdata, sdesign, slogFC, sFDR, NULL, FALSE)

    # Write results out as excel file
    if (!is.null(res.sDE)) {
      #sDE.out <- ExcelOutput(res.sDE, paste0(prefix, out.name))
      sDE.out <- TextOutput(res.sDE, paste0(prefix,out.name))
      rownames(sDE.out) <- NULL
      res.sDE.names <- unique(sDE.out$name)
    } else {
      cat("No signficant DE/DA hits found for analysis of second dataset. Check output file from differential expression analysis. Check your cut-off for differential expression analysis, it may be that these are too stringent.")
    }
    setwd("..")
  }




  if (!is.null(sdata) && technology[2] == "seq") {
    cnames <- colnames(sdata$E)
    sdata <- data.frame(sdata$E)
    colnames(sdata) <- cnames
  }


  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  #                                                                                       ## LASSO Regression ###
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  if (!is.null(lasso)) {
    if(lasso <= 0.0 || lasso > 1.0 ) {
      stop("\n- The input for argument lasso denotes hyperparameter alpha. This value must be set to 0.0 < x < 1.0 for Elastic Net (0.5 is default) or to 1.0 for LASSO regression. Re-run the pipeline again with correct lasso input or remove lasso all together.\n")
    }


    # Length of each group
    len <- as.numeric(table(group))
    test.train <- unique(len >= 19)
    too.few <- unique(len < 9)


    # Stop Lasso if too few samples
    if (TRUE %in% too.few) {
      stop("\n- LASSO cannot be performed, too few samples per group, minimum is 10!\n")
    }


    dir.create("LASSOResults")
    setwd("LASSOResults/")

    group.LASSO <- group
    seeds <- sample(1:1000, 10)
    LASSO.res <- list()

    cat("Cross-validation for grouped multinomial LASSO is running with 10 random seeds, this will take some minutes...")

    if(FALSE %in% test.train) {
      if (databatch == TRUE) {
        if (length(levels(as.factor(group.LASSO))) > 2) {
          for (idx in 1:length(seeds)) {
            LR <- LASSOFeature(seeds[[idx]], data.batch, group.LASSO, lasso, FALSE ,TRUE)
            LASSO.res[[idx]] <-  LR
          }
        } else {
          for (idx in 1:length(seeds)) {
            LR <- LASSOFeature(seeds[[idx]], data.batch, group.LASSO, lasso, FALSE, FALSE)
            LASSO.res[[idx]] <-  LR
          }
        }
      } else {
        if (length(levels(as.factor(group.LASSO))) > 2) {
          for (idx in 1:length(seeds)) {
            LR <- LASSOFeature(seeds[[idx]], data, group.LASSO, lasso, FALSE ,TRUE)
            LASSO.res[[idx]] <-  LR
          }
        } else {
          for (idx in 1:length(seeds)) {
            LR <- LASSOFeature(seeds[[idx]], data, group.LASSO, lasso, FALSE ,FALSE)
            LASSO.res[[idx]] <-  LR
          }
        }
      }
    } else {
      if (databatch == TRUE) {
        if (length(levels(as.factor(group.LASSO))) > 2) {
          for (idx in 1:length(seeds)) {
            LR <- LASSOFeature(seeds[[idx]], data.batch, group.LASSO, lasso, TRUE ,TRUE)
            LASSO.res[[idx]] <-  LR
          }
        } else {
          for (idx in 1:length(seeds)) {
            LR <- LASSOFeature(seeds[[idx]], data.batch, group.LASSO, lasso, TRUE, FALSE)
            LASSO.res[[idx]] <-  LR
          }
        }
      } else {
        if (length(levels(as.factor(group.LASSO))) > 2) {
          for (idx in 1:length(seeds)) {
            LR <- LASSOFeature(seeds[[idx]], data, group.LASSO, lasso, TRUE ,TRUE)
            LASSO.res[[idx]] <-  LR
          }
        } else {
          for (idx in 1:length(seeds)) {
            LR <- LASSOFeature(seeds[[idx]], data, group.LASSO, lasso, TRUE ,FALSE)
            LASSO.res[[idx]] <-  LR
          }
        }
      }
    }


    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Extract results of 10 runs - Write out and plot results
    VarsSelect <- Reduce(intersect, lapply(LASSO.res, '[[', 1))

    if (length(VarsSelect) < 2) {
      stop("\n- There is no overlap in 10 elastic net runs. If you ran LASSO (lasso was et to 1.0) you can try and relax alpha and perform elastic net instead (0.0 < lasso < 1.0). Otherwise you data may have to high of a noise ratio to sample size, LASSO should not be performed.\n")
    }


    VarsSelect <- data.frame(VarsSelect[-1])
    colnames(VarsSelect) <- c("LASSO.Var.Select")


    # Write out LASSO/EN results
    write.table(VarsSelect, paste0(prefix,"_LASSO.txt"), sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Consensus DEA and LASSO
    consensus <- DE.out[DE.out$name %in% VarsSelect$LASSO.Var.Select,]

    if (nrow(consensus) > 0) {
      write.table(consensus, paste0(prefix,"_DEA_LASSO_Consensus.txt"), sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
      pdf(paste0(prefix, "_overlap_DEAA_LASSO_EN.pdf"), height=8, width=12)
      if (length(levels(group)) == 2) {
        venn <- venn.diagram(list(A=unique(as.character(DE.out[DE.out$dir =="up",]$name)), B=unique(as.character(DE.out[DE.out$dir =="down",]$name)), C=as.character(VarsSelect$LASSO.Var.Select)), category.names = c("D(EA) Analysis Up", "D(EA) Analysis Down", "LASSO/EN Regression"), filename=NULL, lwd = 0.7, cat.pos=0, sub.cex = 2, cat.cex= 1.5, cex=1.5, fill=viridis(3, begin = 0.2, end = 0.8, option="cividis"))

      } else {
        venn <- venn.diagram(list(A=unique(as.character(DE.out$name)), B=as.character(VarsSelect$LASSO.Var.Select)), category.names = c("D(EA) Analysis All", "LASSO/EN Regression"), filename=NULL, lwd = 0.7, cat.pos=0, sub.cex = 2, cat.cex= 1.5, cex=1.5, fill=viridis(2, begin = 0.2, end = 0.8, option="cividis"))

      }
      grid.draw(venn)
      dev.off()
    } else {
      cat("There is no consensus between LASSO regression and DEA/DAA.")
    }



    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Cross Validation errors
    LassoRun <- paste0(rep("Run", 10), 1:10)
    CrossValErrormean <- round(unlist(lapply(LASSO.res, '[[', 2)), digits = 4)
    cat(paste0("\nThe average leave-one-out cross validation error for LASSO/elastic-net was: ", mean(CrossValErrormean), "% and the higest error returned from any of the 10 runs was: ", max(CrossValErrormean),"%. Generally the cross validation error should be low ~ 5.0 %, as large errors may indicate a poor model and/or very heterogeneous data. On the other hand, an error of 0 might indicate over-fitting. See CAMPP manual for specifics.\n\n"))
    pCVEM <- data.frame(cbind(CrossValErrormean, LassoRun))
    pCVEM <- ggplot(data=pCVEM, aes(x=LassoRun, y=CrossValErrormean)) + geom_bar(aes(fill = as.factor(LassoRun)), stat="identity") + theme_minimal() + scale_x_discrete(limits=c(LassoRun)) + scale_fill_viridis(begin = 0.0, end = 0.0, discrete=TRUE, option="cividis" ) + theme(legend.position="none") + ylab("CrossValErrormean in %") + theme(axis.text = element_text(size=14), axis.title = element_text(size=16))
    ggsave(paste0(prefix, "_CrossValidationPlot.pdf"), plot = pCVEM)




    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Area under the curve AUC
    if(TRUE %in% test.train) {

      ll <- list()
      llev <- levels(as.factor(group.LASSO))

      for (idx in 1:length(llev)) {
        pos <- which(group.LASSO == as.character(llev[idx]))
        ll[[idx]] <- pos
      }

      my.samp <- unlist(lapply(ll, function(x) sample(x, ceiling((length(x)/4)))))


      if (databatch == TRUE) {
        LASSO.data <- data.batch
      } else {
        LASSO.data <- data
      }


      testD <- data.frame(t(LASSO.data[rownames(LASSO.data) %in% as.character(VarsSelect$LASSO.Var.Select), my.samp]))
      testG <- as.integer(group.LASSO[my.samp])

      trainD <- data.frame(t(LASSO.data[rownames(LASSO.data) %in% as.character(VarsSelect$LASSO.Var.Select), -my.samp]))
      trainG <- as.integer(group.LASSO[-my.samp])

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



  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ### Plotting Results Heatmap ###
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  if (!is.null(plotheatmap)) {
    # color scheme
    colors.hm <- GetColors(as.character(group), colors)

    if (databatch == TRUE){
      hm <- data.batch
      name <- "_batchcorr"
    } else {
      hm <- data
      name <- ""
    }


    if (plotheatmap == "ALL") {
      stop("\nOption ALL is not allowed for heatmap, too heavy! Pick either 'DE', 'DA', 'LASSO', 'EN' or 'Consensus'.\n")
    } else if (plotheatmap %in% c("DA", "DE")) {
      hm <- hm[rownames(hm) %in% res.DE.names,]
    } else {
      if(!is.null(lasso)) {
        if (plotheatmap == "Consensus") {
          hm <- hm[rownames(hm) %in% as.character(consensus$name),]
        }
        if (plotheatmap %in% c("EN", "LASSO")) {
          hm <- hm[rownames(hm) %in% as.character(VarsSelect$LASSO.Var.Select),]
        }
      } else {
        stop("You have specified argument plotheatmap which will produce a heatmap with results from either DA/DE analysis, LASSO/Elastic-Net Regression or the Consensus of these. Input to argument plotheatmap must be a string specifying which results to plot, options are: DA, DE, LASSO, EN, Consensus")
      }
    }


    # heatmap colors in blue
    hm.gradient <- viridis(300, option="cividis")
    range <- c(round(min(DE.out$logFC)), round(max(DE.out$logFC)))

    # Heatmap as pdf
    MakeHeatmap(hm, hm.gradient, colors.hm, colors, groups, paste0(prefix,name), range)

    rm(hm)

  } else {
    cat("\n- No heatmap requested.\n")
  }





  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------





  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ### CORRELATION ANALYSIS ###
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




  # Data (TIF) and sdata correlations

  if (!is.null(sdata) & !is.null(corrby)) {

    dir.create("CorrelationResults")
    setwd("CorrelationResults/")


    retainedSamp <- intersect(colnames(data), colnames(sdata))

    if (corrby == "ALL") {
      retainedvar <- intersect(rownames(data), rownames(sdata))
    } else if (corrby %in% c("DA", "DE")) {
      retainedvar <- intersect(res.DE.names, rownames(sdata))
    } else {
      if(!is.null(lasso)) {
        if (survival == "Consensus" ) {
          retainedvar <- intersect(as.character(consensus$name), rownames(sdata))
        } else {
          retainedvar <- intersect(as.character(VarsSelect$LASSO.Var.Select), rownames(sdata))
        }
      } else {
        stop("You have chosen to perform correlation analysis with results from LASSO/EN BUT this analysis has not been performed. Please re-run pipeline with parameter lasso (see Manual or help (-h))")
      }
    }


    if (databatch == TRUE) {
      data.corr <- data.batch[rownames(data.batch) %in% retainedvar, colnames(data.batch) %in% retainedSamp]
    } else {
      data.corr <- data[rownames(data) %in% retainedvar, colnames(data) %in% retainedSamp]
    }

    if (sdatabatch == TRUE) {
      corr.out <- sdata.batch[rownames(sdata.batch) %in% retainedvar, colnames(sdata.batch) %in% retainedSamp]
    } else {
      corr.out <- sdata[rownames(sdata) %in% retainedvar, colnames(sdata) %in% retainedSamp]
    }

    # Perform correlation analysis and generate overall correlation plot
    res.corr <- CorrAnalysis(data.corr, corr.out, prefix)


    # print out significant hits in Excel
    res.corr$sig.corr <- ifelse(res.corr$fdr <= FDR, "yes", "no")

    write.table(res.corr, paste0(prefix,"_corr.txt"), sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


    # Individual correlation plots
    corr.features <- as.character(res.corr[which(res.corr$sig.corr == "yes"),]$name)
    if (length(corr.features) > 1) {
      CorrelationPlots(data.corr, corr.out, corr.features, prefix)
    } else {
      cat("\n- No significant correlations, no plotting.\n")
    }
    setwd("..")
    try(rm(res.corr,corr.out,data.corr), silent=T)
  }





  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------






  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  ### SURVIVAL ANALYSIS ###
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------






  ### SETTING UP DATA ###
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  if (!is.null(survival)) {

    if (length(intersect(colnames(metadata), must.contain)) < 3) {
      stop("Survival analysis requested, but one or more columns; 'survival', 'outcome', 'outcome.time' are missing from metadata file!")
    }


    dir.create("SurvivalResults")
    setwd("SurvivalResults/")

    # Setting up dataset for survival analysis, matching samples either batch or no batch. Only DA/DE featrues are used.
    metadata.surv <- metadata[which(metadata$survival == 1),]
    samples <- metadata.surv$ids

    if (databatch == TRUE) {
      data.surv <- data.batch[,colnames(data.batch) %in% samples]
    } else {
      data.surv <- data[,colnames(data) %in% samples]
    }



    if (survival == "ALL") {
      cat("\nYou have chosen to use all varibles for survival analysis, this is NOT advisable! See options for survival in with the help argument.\n")
    } else if (survival %in% c("DA", "DE")) {
      data.surv <- data.surv[rownames(data.surv) %in% res.DE.names,]
    } else {
      if(!is.null(lasso)) {
        if (survival == "Consensus" ) {
          data.surv <- data.surv[rownames(data.surv) %in% as.character(consensus$name),]
        } else {
          data.surv <- data.surv[rownames(data.surv) %in% as.character(VarsSelect$LASSO.Var.Select),]
        }
      } else {
        stop("You have chosen to perform survival analysis with results from LASSO/EN but this analysis has not been performed. Please re-run pipeline with parameter lasso (see Manual)")
      }
    }





    # If user has specified covariates for survival analysis, extract these.
    if(is.null(covarS)) {
      surv_object <- data.frame(t(data.surv), as.numeric(metadata.surv$outcome.time), metadata.surv$outcome)
      colnames(surv_object) <- c(rownames(data.surv), "outcome.time", "outcome")
    } else {
      my.covars <- metadata.surv[,colnames(metadata.surv) %in% covarS]
      surv_object <- data.frame(t(data.surv), as.numeric(metadata.surv$outcome.time), metadata.surv$outcome, my.covars)
      colnames(surv_object) <- c(rownames(data.surv), "outcome.time", "outcome", covarS)
    }


    # datadist format for cph function with cubic splines.
    features <- rownames(data.surv)
    dd <- datadist(surv_object)
    options(datadist=dd)

    mylist <- list(features, dd, surv_object)




    ### User Specified Covariates ###
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




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
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------




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


    # # Colnames with and without user-specified covariates.
    # if(is.null(covarS)) {
    #     colnames(covariate_linearity) <- c("Feature")
    # } else {
    #     covar.names <- unlist(strsplit(covarS, split="[+]"))
    #     if(length(grep("rcs", covar.names)) > 0) {
    #         covar.names <- c("Feature", covar.names[grep("rcs", covar.names)], "Total")
    #         covar.names <- gsub("rcs\\(|\\)", "", covar.names)
    #         colnames(covariate_linearity) <- covar.names
    #     } else {
    #         colnames(covariate_linearity) <- c("Feature")
    #     }
    # }



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
    # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------



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
      features.surv <- features[-c((length(features)-length(covar)):length(features))]
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
    survival.data <- SurvivalCOX(survival.results, prefix, survplot)

    write.table(survival.data, paste0(prefix,"_survival.txt"), sep = "\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

    setwd("..")
    rm(dd, pha_check, pha.fail.test, survival.results, features.surv, survival.data)
  }






  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # SOURCE FUNCTIONS FROM THE FUNCTIONS SCRIPT - PART 3
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  # rm(DAFeature, DAFeatureApply, LASSOFeature, TextOutput, MakeHeatmap, SurvivalCOX, CorrAnalysis, CorrelationPlots)
  # gc(full = TRUE)
  #


  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  #                                                               ## Weighed Gene Co-Expression Network Analysis ###
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


  if (!is.null(WGCNA)) {

    dir.create("WGCNAResults")
    setwd("WGCNAResults/")

    if (databatch == TRUE){
      data.WGCNA <- t(data.batch)
    } else {
      data.WGCNA <- t(data)
    }

    if (WGCNA %in% c("DA", "DE")) {
      data.WGCNA <- data.WGCNA[,colnames(data.WGCNA) %in% res.DE.names]
    }

    dir.create("WGCNAPlots")
    setwd("WGCNAPlots/")
    # Check data
    gsg <- goodSamplesGenes(data.WGCNA)
    cat(paste0("- Data set is OK for WGCNA - ", gsg$allOK,".\n"))
    WGCNAres <- WGCNAAnalysis(data.WGCNA, cutoffWGCNA, prefix)
    setwd("..")

    if (WGCNA %in% c("DA", "DE")) {
      logFCs <- DE.out
      logFCs$gene <- DE.out$name
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





  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------






  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  #                                                               ## Protein-Protein and MiRNA-Gene Network Analysis ###
  # ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------







  if (!is.null(PPI)) {

    DB <- DownloadPPInt(PPI[[1]])

    dir.create("InteractionResults")
    setwd("InteractionResults/")
    dir.create("InteractionPlots")

    DE.out$name <- gsub("_", "-", DE.out$name)

    PPIs <- GetDeaPPInt(DB, DE.out)


    if (!is.null(GmiRI)) {

      sDE.out$name <- gsub("_", "-", sDE.out$name)

      cat("\nSearching for miRNA-gene interaction - this may take up to 10 min!\n")
      GmiRIs <- GetGeneMiRNAInt(GmiRI[[2]], sDE.out, DE.out)

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

    DE.out$name <- gsub("_", "-", DE.out$name)

    cat("\nSearching for miRNA-gene interaction - this may take up to 10 min!\n")
    GmiRIs <- GetGeneMiRNAInt(GmiRI[[2]], DE.out)
    GmiRTrim <- TrimWriteInt(GmiRIs)

    setwd("InteractionPlots/")
    cat("\nPlotting set-wise networks - this may take up to 10 min!\n")
    PlotInt(GmiRTrim)

    setwd("../..")

  } else {
    cat("\nNo Interaction Networks requested\n.")
  }

  setwd("../")
  cat("\nCAMPP RUN DONE!\n")

}

