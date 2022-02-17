library(CAMPP2)


getwd()
setwd("/Users/nikto/opt/campp_bioconducor/CAMPP2/data/BRCA")



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

#runCampp2(prefix="results_test1", data=dataset1, metadata=metadata1, groups=c("IDs", "diagnosis"), technology=c("seq"))

# runCampp2(PPint=c("ensembl_gene_id","stringdatabase"),WGCNA="DE", signif=c(1,1),survival="DE",corr=NULL,plotheatmap="DE",batches="tumor_stage",kmeans=TRUE, plotmds=TRUE,prefix="test_dual_surv_covar", data=dataset1, metadata=metadata1, groups=c("IDs", "diagnosis"), technology=c("seq"))
# runCampp2(WGCNA="DE", signif=c(1,1),survival="DE",corr=NULL,plotheatmap="DE",batches="tumor_stage",kmeans=TRUE, plotmds=TRUE,prefix="test_dual_surv_covar2", data=dataset1, metadata=metadata1, groups=c("IDs", "diagnosis"), technology=c("seq"))
# runCampp2(cutoffWGCNA=c(10,25,25),survplot=30, WGCNA="DE", signif=c(1,1),survival="DE",corr=NULL,plotheatmap="DE",batches="tumor_stage",kmeans=TRUE, plotmds=TRUE,prefix="test_dual_surv_covar3", data=dataset1, metadata=metadata1, groups=c("IDs", "diagnosis"), technology=c("seq"))
# runCampp2(covar=c(FALSE,"age"),cutoffWGCNA=c(10,25,25),survplot=30, WGCNA="DE", signif=c(1,1),survival="DE",corr=NULL,plotheatmap="DE",batches="tumor_stage",kmeans=TRUE, plotmds=TRUE,prefix="test_dual_surv_covar4", data=dataset1, metadata=metadata1, groups=c("IDs", "diagnosis"), technology=c("seq"))
# runCampp2(PPint=c("ensembl_gene_id","stringdatabase"), prefix="test_dual_surv_covar5", data=dataset1, metadata=metadata1, groups=c("IDs", "diagnosis"), technology=c("seq"))
# runCampp2(lasso=0.5, prefix="test_dual_surv_covar6", data=dataset1, metadata=metadata1, groups=c("IDs", "diagnosis"), technology=c("seq"))
#
# runCampp2(WGCNA="DE", signif=c(1,1,1,1),survival="DE",corr=NULL,plotheatmap="DE",batches=c("tumor_stage","tumor_stage"),kmeans=TRUE, plotmds=TRUE,prefix="test_dual_surv_covar7", data=dataset1, sdata=dataset2, metadata=metadata1,smetadata=metadata2, groups=c("IDs", "diagnosis","IDs", "diagnosis"), technology=c("seq","seq"))

data=dataset1
metadata=metadata1
groups=c("IDs", "diagnosis")
technology=c("seq")
sdata<-"sdata"
prefix<-"prefix_test"



runCampp2_dev(data=data, metadata=metadata, groups=groups, technology=technology, prefix=prefix)




