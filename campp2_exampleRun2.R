library(CAMPP2)
library(edgeR)
library(sva)
library(pROC)
library(VennDiagram)

getwd()
setwd("/data/user/mathilde/CAMPP2/data_glycans")


data_glycan1<-as.data.frame(importCounts("glycandata.xlsx"))
data_Sglycan1<-as.data.frame(importCounts("glycanSdata.xlsx"))

dataset1<-merge(x=data_glycan1,y=data_Sglycan1,by='row.names', all=TRUE)
rownames(dataset1)<-dataset1[,1]
dataset1<-dataset1[,-1]
metadata$IDs<-rownames(metadata)

#prepare metadata
metadata<-as.data.frame(importCounts("glycanmetadata.xlsx"))

#Check sample names in data meet sample names in metadata; The sample names must match
head(colnames(dataset1))
head(rownames(metadata))

#if needed, correct samples names
keep_samples1<-gsub("\\.","-",colnames(dataset1))
colnames(dataset1)<-keep_samples1

#subset metadata only to samples in included in a appropriate dataset
metadata1<-subset(metadata, metadata$IDs %in% keep_samples1)

#sort dataset1 and metadata by ID - THIS IS ABSOLUTELY NECESSARY
dataset1<-dataset1[ , order(names(dataset1))]
metadata1<-metadata1[order(metadata1$IDs),]

setwd("../")
getwd()

#runCampp2(WGCNA="DE", signif=c(1,1,1,1),survival="DE",plot.heatmap="DE",batches=c("tumor_stage","tumor_stage"),kmeans=TRUE, plot.mds=TRUE,prefix="test_dual_surv_covar7", data1=dataset1, data2=dataset2, metadata1=metadata1,metadata2=metadata2, groups=c("IDs", "diagnosis","IDs", "diagnosis"), technology=c("seq","seq"),group1)
runCampp2(signif=c(1,1,1,1),survival="DE",correlation=FALSE,plot.heatmap="DE",batches=c("tumor_stage","tumor_stage"), plot.mds=TRUE, data1=dataset1, data2=NULL, metadata1=metadata1,metadata2=NULL, groups=c("IDs","Diagnosis"), technology=c("seq","seq"))

