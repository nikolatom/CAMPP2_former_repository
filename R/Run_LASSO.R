#' @title A wrapper function for LASSO/Elastic net/Ridge regression
#' @description The function runs LASSO/Elastic net/Ridge regression (depending
#' on a hyperparameter alpha value. For Ridge regression, alpha = 0;
#' for Elastic net, alpha >0 & <1; for LASSO, alpha = 1).
#' The regression is calculated 10 times for 10 random seeds and the intersect
#' of the significant features (e.g., genes) is obtained.
#' Parameter Lambda is automatically
#' estimated using k-fold cross-validation (cv.glmnet) by default.
#' The function optionally (>20 samples) calculates a double
#' cross-validation using a validation set which represents 1/4 of the samples
#' from each group.
#' Results are provided in the form of:
#' 1) the names of the significant features passing filters based on their
#' coefficient values
#' 2) classification error rate (datasets of >20 samples)
#' 3) roc.res - AUC value (datasets of >20 samples)
#' 4) cross validation error plot
#' @param data a data frame of expression/abundance counts. It's recommended
#' to use normalized and batch corrected data.
#' @param group a factor specifying group for each sample (e.g. could be
#' represented by a column from a metadata file)
#' @param alpha a numeric vector specifying hyperparameter alpha. This value
#' must be set to 0.0 < x < 1.0 for Elastic Net (0.5 is default) or to 1.0
#' for LASSO regression or to 0.0 for Ridge regression.
#' @param min.coef a numeric vector specifying a threshold for features'
#' filtering (e.g. genes) based on the coefficients which are calculated during
#' model fitting. Default value is > 0.
#' @param nfolds a numeric vector describing number of folds during Lambda
#' estimation which is based on a cross-validation. Although nfolds can be
#' as large as the sample size (leave-one-out CV), it is not recommended for
#' large datasets. Smallest value allowable is nfolds=3. Default is 10.
#' @param prefix a character string defining a prefix of output file.
#' @export
#' @import glmnet
#' @import parallel
#' @import doMC
#' @import viridis
#' @import pROC
#' @import ggplot2
#' @import grid
#' @import nnet
#' @return a list of:
#' 1) VarSelect - a matrix of feature names having coefficients of best model
#' passing the filters (threshold defined by min.coef)
#' 2) cv.error - logical/numeric value describing miss classification error
#' rate (datasets of >20 samples)
#' 3) roc.res - AUC value (datasets of >20 samples)
#' 4) cross validation error plot
#' @examples {
#' ###LASSO can be run on data where there is at least 8 samples
#' ###per group. Here we create a large dataset
#' campp2_test_data_LASSO<-cbind(campp2_brca_1,campp2_brca_2)
#' campp2_test_metadata_LASSO<-rbind(campp2_brca_1_meta, campp2_brca_2_meta)
#' campp2_test_data_LASSO_replaceNAs<-ReplaceNAs(data=campp2_test_data_LASSO)
#' campp2_test_data_LASSO_zeroFix<-FixZeros(data=campp2_test_data_LASSO_replaceNAs,group=campp2_test_metadata_LASSO$diagnosis, remove.sparse.features=TRUE)
#' campp2_test_data_LASSO_normalized<-NormalizeData(data=campp2_test_data_LASSO_zeroFix,group=campp2_test_metadata_LASSO$diagnosis,standardize="TMM",transform="voom",technology="seq")
#' campp2_test_data_LASSO_batchCorrected<-BatchCorrect(data=campp2_test_data_LASSO_normalized,batch=campp2_test_metadata_LASSO$tumor_stage,group=campp2_test_metadata_LASSO$diagnosis,technology="seq")
###run lasso on a large dataset
#' runLASSO(data=campp2_test_data_LASSO_batchCorrected, group=campp2_test_metadata_LASSO$diagnosis, alpha=0.5, min.coef=0, nfolds=10, prefix="test")
#' }


runLASSO <- function(data, group, alpha, min.coef=0, nfolds=10, prefix=NULL){

    ##Checking parameter alpha
    if(alpha <= 0.0 || alpha > 1.0 ) {
        stop("\n- The input for argument lasso denotes hyperparameter alpha. This value must be set to 0.0 < x < 1.0 for Elastic Net (0.5 is default) or to 1.0 for LASSO regression. Re-run the pipeline again with correct lasso input or remove lasso all together.\n")
    }else if (alpha == 0.0 ){
        print("According to parameter alpha, Ridge regression will be run.")
    }else if (alpha > 0.0 && alpha < 1){
        print("According to parameter alpha, Elastic net will be run.")
    }else if (alpha == 1){
        print("According to parameter alpha, LASSO will be run.")
    }


    # Length of each group
    group.size <- as.numeric(table(group))
    test.train <- unique(group.size >= 20)
    if(FALSE %in% test.train){
        print("Double cross validation is not going to be done due to the low number (<20) of samples in at least one of the sample groups.")
    }


    # Stop Lasso if too few samples
    too.few <- unique(group.size <= 7)
    if (TRUE %in% too.few) {
        stop("\n- LASSO cannot be performed, too few samples per group, minimum is 8!\n")
    }


    ###Setting up the seeds
    seeds <- sample(1:1000, 10) #generate 10 random seeds between 1:1000
    LASSO.res <- list()

    cat("Cross-validation for grouped multinomial LASSO/EN/Ridge regression is running with 10 random seeds.")

    ### test.train could be TRUE FALSE if one group is small 8-19; validation would not be done (validation sets not created)
    if(FALSE %in% test.train) {
        for (idx in 1:length(seeds)) {
            LR <- LASSOFeature(seed=seeds[[idx]], data, group=group,alpha, validation=FALSE, min.coef=0, nfolds = 10)
            LASSO.res[[idx]] <-  LR
        }
    } else {   ###with validation (>=20 samples in each group)
        for (idx in 1:length(seeds)) {
            LR <- LASSOFeature(seed=seeds[[idx]], data, group=group,alpha, validation=TRUE, min.coef=0, nfolds = 10)
            LASSO.res[[idx]] <-  LR
         }
    }


    ### Extract results from 10 runs and
    VarsSelect <- Reduce(intersect, lapply(LASSO.res, '[[', 1)) ##Creates an intersect of significant features from 10 data frames saved in LASSO.res (results are different because of different seeds)

    if (length(VarsSelect) < 2) {
        stop("\n- There is no overlap in 10 LASSO/elastic net/Ridge regression runs. If you ran LASSO (lasso was set to 1.0) you can try and relax alpha and perform elastic net instead (0.0 < lasso < 1.0). Otherwise you data may have to high of a noise ratio to sample size, LASSO should not be performed.\n")
    }


    VarsSelect <- data.frame(VarsSelect[-1])
    colnames(VarsSelect) <- c("VarsSelect")


    ### Cross Validation errors
    if(FALSE %in% test.train){
        cv.error=NA
        print("Cross-validation error rate using test set wasn't calculated due to the small size of the sample group(s)")
    }else{
        run.number <- as.character(1:10)
        cross.val.error.mean <- round(unlist(lapply(LASSO.res, '[[', 2)), digits = 4) ##print error means from the results
        cat(paste0("The average (the analysis was run with 10 random seeds) cross validation error for LASSO/Elastic net/Ridge regression was: ", mean(cross.val.error.mean), "% and the higest error returned from any of the 10 runs was: ", max(cross.val.error.mean),"%. Generally the cross validation error should be low ~ 5.0 %, as large errors may indicate a poor model and/or very heterogeneous data. On the other hand, an error of 0 might indicate over-fitting. \n\n"))
        ## originally, it was saying "leave-one-out cross validation error" - I think it's not true
        cv.error <- data.frame(cbind(cross.val.error.mean, run.number)) #cv.error - percentage of cross validation error
        cv.error.plot <- ggplot(data=cv.error, aes(x=run.number, y=cross.val.error.mean)) + geom_bar(aes(fill = as.factor(run.number)), stat="identity") + theme_minimal() + scale_x_discrete(limits=c(run.number)) + scale_fill_viridis(begin = 0.0, end = 0.0, discrete=TRUE, option="cividis" ) + theme(legend.position="none") + ylab("cross.val.error.mean in %") + theme(axis.text = element_text(size=14), axis.title = element_text(size=16))
        ggsave(paste0(prefix, "_CrossValidationPlot.pdf"), plot = cv.error.plot)

    }


    ### Area under the curve AUC performed only if all the sample groups are large enough
    if(!FALSE %in% test.train) {

        ll <- list()
        llev <- levels(as.factor(group))
        for (idx in 1:length(llev)) {
            pos <- which(group == as.character(llev[idx]))
            ll[[idx]] <- pos
        }

        samp <- unlist(lapply(ll, function(x) sample(x, ceiling((length(x)/4)))))

        testD <- data.frame(t(data[rownames(data) %in% as.character(VarsSelect[,1]), samp]))
        testG <- as.character(group[samp])

        trainD <- data.frame(t(data[rownames(data) %in% as.character(VarsSelect[,1]), -samp]))
        trainG <- as.character(group[-samp])

        mn.net <- nnet::multinom(trainG ~ ., data=trainD)  ##is multinom correct?
        mn.pred <- predict(mn.net, newdata=testD, type="prob")
        roc.res <- multiclass.roc(testG, mn.pred)
        roc.res <- data.frame(round(as.numeric(sub(".*: ", "", roc.res$auc)), digits = 2))

        colnames(roc.res) <- "AUC"

        cat(paste0("The area under the curve (AUC) for variables selected (based on intersection) from 10 LASSO/EN/Rigge regression runs was: ", roc.res,". \n"))
    } else{
        roc.res=NA
    }

    return(list("VarsSelect"=VarsSelect, "cv.error"=cv.error, "roc.res"=roc.res))

}


