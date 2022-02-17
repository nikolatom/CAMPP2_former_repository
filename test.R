library(stringr)
library(dplyr)

signif=c(1,1)
data=tumor_testData1
kmeans=NULL
metadata=metadata1
technology=c("seq")
groups=c("IDs", "diagnosis")


sdata=dataset2
smetadata=metadata2

batches=c("tumor_stage","tumor_stage")
datacheck=NULL
standardize=NULL
transform=NULL
plotmds=NULL
plotheatmap=NULL
kmeans=NULL
colors=NULL
prefix="test_surv"
corr=NULL
lasso=NULL
WGCNA=NULL
cutoffWGCNA=NULL
survival="DE"
covar=NULL
stratify=NULL
survplot=NULL
PPint=NULL
GenemiRint=NULL

samples_test<-colnames(tumor_testData1)
metadata_test<-as.data.frame(metadata[metadata$IDs %in% samples_test,])

colnames(tumor_testData1)<-gsub("-",".",samples_test)
metadata_test$IDs<-as.data.frame(gsub("-",".",metadata_test$IDs))
colnames(metadata_test[,1])<-"IDs"

test_run <- function (data, metadata, sdata=NULL, smetadata=NULL, technology, groups, batches=NULL, datacheck=TRUE, standardize=NULL, transform=NULL, plotmds=NULL, plotheatmap=NULL, kmeans=NULL, signif=NULL, colors=NULL, prefix=NULL, corr=NULL, lasso=NULL, WGCNA=NULL, cutoffWGCNA=NULL, survival=NULL, covar=NULL, stratify=NULL, survplot=NULL, PPint=NULL, GenemiRint=NULL){
    parseArguments(data, metadata, sdata, smetadata, technology, groups, batches, datacheck, standardize, transform, plotmds, plotheatmap, kmeans, signif, colors, prefix, corr, lasso, WGCNA, cutoffWGCNA, survival, covar, stratify, survplot, PPint, GenemiRint)
}

test_run(data,metadata,technology,groups)
test_run(data=tumor_testData1,metadata=metadata_test,technology=technology,groups=groups)

runCampp2(data=tumor_testData1,metadata=metadata_test,technology="seq",groups=groups)

