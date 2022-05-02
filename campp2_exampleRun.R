library(CAMPP2)
library(edgeR)
library(sva)
library(pROC)
library(VennDiagram)
library(rms)
library(dynamicTreeCut)

getwd()
setwd("../data_nik")



data_normal1<-as.data.frame(importCounts("normal_testData1.txt"))
data_normal2<-as.data.frame(importCounts("normal_testData2.txt"))
data_tumor1<-as.data.frame(importCounts("tumor_testData1.txt"))
data_tumor2<-as.data.frame(importCounts("tumor_testData2.txt"))

dataset1<-merge(x=data_normal1,y=data_tumor1,by='row.names', all=TRUE)
dataset2<-merge(x=data_normal2,y=data_tumor2,by='row.names', all=TRUE)
rownames(dataset1)<-dataset1[,1]
dataset1<-dataset1[,-1]
rownames(dataset2)<-dataset2[,1]
dataset2<-dataset2[,-1]

#prepare metadata using getherMetadata.R
source('./gatherMetadata.R')
#remove spaces from metadata columns
metadata$diagnosis<-gsub(" ","_",metadata$diagnosis)
metadata$tumor_stage<-gsub(" ","_",metadata$tumor_stage)

#Check sample names in data meet sample names in metadata; The sample names must match
head(colnames(dataset1[-1]))
head(metadata$IDs)

#if needed, correct samples names
keep_samples1<-gsub("\\.","-",colnames(dataset1))
keep_samples2<-gsub("\\.","-",colnames(dataset2))
colnames(dataset1)<-keep_samples1
colnames(dataset2)<-keep_samples2

#subset metadata only to samples in included in a appropriate dataset
metadata1<-subset(metadata, metadata$IDs %in% keep_samples1)
metadata2<-subset(metadata, metadata$IDs %in% keep_samples2)

#sort dataset1 and metadata by ID - THIS IS ABSOLUTELY NECESSARY
dataset1<-dataset1[ , order(names(dataset1))]
metadata1<-metadata1[order(metadata1$IDs),]
dataset2<-dataset2[ , order(names(dataset2))]
metadata2<-metadata2[order(metadata2$IDs),]
dataset3<-cbind(dataset1,dataset2)
metadata3<-rbind(metadata1,metadata2)

runCampp2(prefix="test_normalization", data1=dataset1, data2=dataset2, metadata1=metadata1,metadata2=metadata2, groups=c("IDs", "diagnosis","IDs", "diagnosis"), technology=c("seq","seq"))
