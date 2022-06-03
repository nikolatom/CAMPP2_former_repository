library(CAMPP2)
library(edgeR)
library(sva)
library(pROC)
library(VennDiagram)
library(rms)
library(dynamicTreeCut)

###This command will run CAMPP2 with default settigs using example data which are included in the package
runCampp2(batches=c("tumor_stage","tumor_stage"),prefix="test_CAMPP2", data1=dataset1, data2=dataset2, metadata1=metadata1,metadata2=metadata2, groups=c("IDs", "diagnosis","IDs", "diagnosis"), technology=c("seq","seq"))

###Code for testing the function for Replacing NAs
n<-1000
Na_dataset1<-apply (dataset1, 2, function(x) {x[sample( c(1:n), floor(n/10))] <- NA; x} )
ReplaceNAs(Na_dataset1)

###Code for checking zeros will be added when finished
###Code for normalization (finished) will be added when function for checking the zeros will be on place
