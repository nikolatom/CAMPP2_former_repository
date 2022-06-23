library(CAMPP2)
library(edgeR)
library(sva)
library(pROC)
library(VennDiagram)
library(rms)
library(dynamicTreeCut)

###This command will run CAMPP2 with default settigs using example data which are included in the package
runCampp2(batches=c("tumor_stage","tumor_stage"),prefix="test_CAMPP2", data1=dataset1, data2=dataset2, metadata1=metadata1,metadata2=metadata2, groups=c("IDs", "diagnosis","IDs", "diagnosis"), technology=c("seq","seq"))

###An example code for testing the function for Replacing NAs
n<-10000
Na_dataset1<-apply (dataset1, 2, function(x) {x[sample( c(1:n), floor(n/10))] <- NA; x} )
replacedNAs<-ReplaceNAs(data=Na_dataset1)

###An example code for fixing zeros and negative values
zerofix<-FixZeros(data=replacedNAs,group=metadata1$diagnosis,data, group, remove.sparse.features=TRUE)

###An example code for normalization
normalized_data<-NormalizeData(data=zerofix,group=metadata1$diagnosis,standardize="none",transform="none",technology="seq")

###An example code for batch correction
batch_corrected_data<-BatchCorrect(normalized_data,metadata1$tumor_stage,metadata1$diagnosis,technology="seq")