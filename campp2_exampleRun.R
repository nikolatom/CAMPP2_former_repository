library(CAMPP2)
library(edgeR)
library(sva)
library(pROC)
library(VennDiagram)
library(rms)
library(dynamicTreeCut)

###This command will run CAMPP2 with default settigs using example data which are included in the package
runCampp2(batches=c("tumor_stage","tumor_stage"),prefix="test_CAMPP2", data1=campp2_brca_1, data2=campp2_brca_2, metadata1=campp2_brca_1_meta, metadata2=campp2_brca_2_meta, groups=c("IDs", "diagnosis","IDs", "diagnosis"), technology=c("seq","seq"))

###An example code for testing the function for Replacing NAs
n<-10000
Na_campp2_brca_1<-apply (campp2_brca_1, 2, function(x) {x[sample( c(1:n), floor(n/10))] <- NA; x} )
replacedNAs<-ReplaceNAs(data=Na_campp2_brca_1)

###An example code for fixing zeros and negative values
zerofix<-FixZeros(data=replacedNAs,group=campp2_brca_1_meta$diagnosis, remove.sparse.features=TRUE)

###An example code for normalization
normalized_data<-NormalizeData(data=zerofix,group=campp2_brca_1_meta$diagnosis,standardize="none",transform="none",technology="seq")

###An example code for batch correction
batch_corrected_data<-BatchCorrect(normalized_data,campp2_brca_1_meta$tumor_stage,campp2_brca_1_meta$diagnosis,technology="seq")

###An example code for running differential gene expression analysis
runCampp2(prefix="Testing_DEA", data1=campp2_brca_1, metadata1=campp2_brca_1_meta, groups=c("IDs", "diagnosis"), technology=c("seq"), block=c(campp2_brca_1_meta$tumor_stage))

###An example code for running differential expression analysis with subtype analysis and visualizations
runCampp2(prefix="Testing_DEA", data1=campp2_brca_1, metadata1=campp2_brca_1_meta, groups=c("IDs", "subtype"), technology=c("seq"), plot.DEA=TRUE, plot.heatmap="DEA")
