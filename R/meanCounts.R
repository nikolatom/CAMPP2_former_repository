#' @title Get mean and SD of the feature counts
#' @description A function for obtaining an average feature counts and SD for
#' each gene across all the samples and for each sample group. Results are
#' provided in a form of a list of data frames and figures which are exported
#' into .pdf.
#' @param data a data frame of feature counts.
#' @param group a factor specifying group for each sample (e.g. could be represented by a column from a metadata file)
#' @export
#' @import ggplot2
#' @return a list of data frames (for each sample group) describing mean and SD
#' for all sample groups and also for each group of samples. Results are
#' visualized via graphs saved into current directory as pdf files.
#' @examples {
#' meanCounts(campp2_brca_1_batchCorrected, campp2_brca_1_meta$diagnosis)
#' }


meanCounts<-function(data, group){

    mean.all<-as.data.frame(sort(rowMeans(data)))
    data.mean.all<-merge(mean.all,data,by=0, all=TRUE, sort=FALSE) #keep sorting from mean.all
    sd.all<-as.data.frame(apply(data.mean.all[,3:ncol(data.mean.all)], 1, sd)) #calculate sd
    rownames(sd.all)<-data.mean.all[,1]
    index.all<-as.data.frame(1:nrow(mean.all)) #create index for each gene
    rownames(index.all)<-data.mean.all[,1]
    mean.sd.all<-merge(mean.all, sd.all, by=0, all=TRUE, sort=FALSE) #merge index, mean and sd based on mean
    rownames(mean.sd.all)<-mean.sd.all[,1]
    mean.sd.all<-mean.sd.all[,-1]
    mean.sd.index.all<-merge(mean.sd.all, index.all, by=0, all=TRUE, sort=FALSE)  ###index, mean and sd for each gene; sorted by mean
    colnames(mean.sd.index.all)<-c("gene_id", "mean.all", "sd.all", "index.all")

    group.level<-levels(as.factor(group))
    group.data<-list()

    for(i in group.level){
        sample.idx.group<-(which(group == i))
        data.group<-data[,sample.idx.group] ##extract samples of interest

        ##############
        mean.group<-as.data.frame(sort(rowMeans(data.group)))
        data.mean.group<-merge(mean.group,data.group,by=0, all=TRUE, sort=FALSE)
        sd.group<-as.data.frame(apply(data.mean.group[,3:ncol(data.mean.group)], 1, sd)) #calculate sd
        rownames(sd.group)<-data.mean.group[,1]
        index.group<-as.data.frame(1:nrow(mean.group)) #create index for each gene
        rownames(index.group)<-data.mean.group[,1]
        mean.sd.group<-merge(mean.group, sd.group,by=0, all=TRUE, sort=FALSE) #merge index, mean and sd based on mean
        rownames(mean.sd.group)<-mean.sd.group[,1]
        mean.sd.group<-mean.sd.group[,-1]
        mean.sd.index.group<-merge(mean.sd.group,index.group,by=0, all=TRUE, sort=FALSE)  ###index, mean and sd for each gene; sorted by mean
        colnames(mean.sd.index.group)<-c("gene_id.group", "mean.group", "sd.group", "index.group")

        df.group<-merge(mean.sd.index.all, mean.sd.index.group,by=1, all=TRUE, sort=FALSE) ###keep sorting from mean.sd.index.all
        group.data[[i]]<-df.group
        #################

        colors <- c("SD.all" = "yellow",  "SD.group" = "violet", "mean.all" = "black", "mean.group" = "blue")

        res.plot<-
            ggplot(data=df.group, aes(x=index.all, y=mean.all)) +
            geom_errorbar(aes(ymin=mean.all-sd.all, ymax=mean.all+sd.all,color="SD.all"),alpha=0.1)+
            geom_errorbar(aes(ymin=mean.group-sd.group, ymax=mean.group+sd.group, color="SD.group") ,alpha = 0.1)+
            geom_point(aes(x=index.all, y=mean.group,color="mean.group"), size=0.1)+
            geom_point(color="black", size=0.1)+
            scale_color_manual(name="",values = colors,  breaks=c('SD.all', 'SD.group', 'mean.all', 'mean.group'))+
            labs(title = "Plot of the mean and SD of the feature counts",y = "feature counts", x = "feature index (sorted by mean of all samples)", color="Legend")+
            scale_alpha(guide = 'none')


        ggsave(paste0("mean_feature_count_",i,".pdf"), width = 4, height = 4)

    }
    return(group.data)
}







