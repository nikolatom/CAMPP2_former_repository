library(dplyr)
setwd("/Users/nikto/opt/campp_bioconducor/CAMPP2/data/BRCA")
#usage:
    ##source function first
    load("./BRCA_dataPrep.rda")
    normal_tumor_counts<-extract_tumor_normal_samples(dataPrep)
    normal_counts_paired<-as.data.frame(normal_tumor_counts$normal_paired)
    tumor_counts_paired<-as.data.frame(normal_tumor_counts$tumor_paired) #including technical replicates
    normal_counts<-as.data.frame(normal_tumor_counts$normal)
    tumor_counts<-as.data.frame(normal_tumor_counts$tumor)
    normal_duplicates<-normal_tumor_counts$N_dup_names
    tumor_duplicates<-normal_tumor_counts$T_dup_names
    all_normal_names<-normal_tumor_counts$all_normal_names
    all_tumor_names<-normal_tumor_counts$all_tumor_names

    write.table(tumor_counts_paired, file="tumor_counts_paired.txt", sep = "\t")
    write.table(normal_counts_paired, file="normal_counts_paired.txt", sep = "\t")
    write.table(tumor_counts, file="tumor_counts.txt", sep = "\t")
    write.table(normal_counts, file="normal_counts.txt", sep = "\t")
    write.table(tumor_duplicates, file="tumor_duplicates.txt", sep = "\t")
    write.table(normal_duplicates, file="normal_duplicates.txt", sep = "\t")

    ###Creating test data:
    normal_testData1<-normal_counts_paired[,1:50]
    normal_testData2<-normal_counts_paired[,51:100]
    tumor_testData1<-tumor_counts_paired[,1:50]
    tumor_testData2<-tumor_counts_paired[,51:100]
    write.table(normal_testData1, file="normal_testData1.txt", sep = "\t")
    write.table(normal_testData2, file="normal_testData2.txt", sep = "\t")
    write.table(tumor_testData1, file="tumor_testData1.txt", sep = "\t")
    write.table(tumor_testData2, file="tumor_testData2.txt", sep = "\t")



    ##This function returns list "normal_tumor_counts" of:
    ##normal_counts_paired
    ##tumor_counts_paired
    ##normal_counts
    ##tumor_counts
    ##normal_duplicates
    ##tumor_duplicates
    ##all normal names
    ##all tumor names

extract_tumor_normal_samples <- function(gene_count_matrix, normal_sampleRegExp = '....-..-....-1[1-9].-...-....-..', tumor_sampleRegExp = '....-..-....-0[1-9].-...-....-..' ){
    #extract tumor and normal samples based on sample tag regular expression
    #Sample type is defined as: Tumor types range from 01 - 09, normal types from 10 - 19 and control samples from 20 - 29 (coded between 3rd and 4th "-").
    normal<-as.data.frame(gene_count_matrix) %>% select(matches(normal_sampleRegExp))
    tumor<-as.data.frame(gene_count_matrix) %>% select(matches(tumor_sampleRegExp))

    #get only data from the participants having both tumor and normal samples
    #participant ID is coded between 2nd and 3rd "-".
    normal_col<-as.matrix(colnames(normal))
    tumor_col<-as.matrix(colnames(tumor))
    normal_list<-strsplit(normal_col,"-")
    tumor_list<-strsplit(tumor_col,"-")

    normal_names<-NULL
    for (ID in normal_list) {
        parts<-unlist(ID)
        normal_names<-append(normal_names,parts[3])
    }
    tumor_names<-NULL
    for (ID in tumor_list) {
        parts<-unlist(ID)
        tumor_names<-append(tumor_names,parts[3])
    }

    #GET DUPLICATED samples (first occurence of sample is not counted)
    tumor_duplicates<-duplicated(tumor_names)
    T_dup_index<-which(tumor_duplicates,arr.ind=FALSE, useNames = TRUE)
    T_dup_names<-NULL
    for (ID in T_dup_index) {
        T_dup_names<-append(T_dup_names,tumor_names[ID])
    }

    normal_duplicates<-duplicated(normal_names)
    N_dup_index<-which(normal_duplicates,arr.ind=FALSE, useNames = TRUE)
    N_dup_names<-NULL
    for (ID in N_dup_index) {
        N_dup_names<-append(N_dup_names,tumor_names[ID])
    }


    intersect_tumor_normal <- intersect(normal_names,tumor_names)
    normal_paired<-as.data.frame(normal %>% select(contains(intersect_tumor_normal)))
    tumor_paired<-as.data.frame(tumor %>% select(contains(intersect_tumor_normal)))
    output<-list(normal_paired,tumor_paired, normal, tumor, N_dup_names, T_dup_names, normal_names, tumor_names)
    names(output) <- c("normal_paired","tumor_paired", "normal", "tumor", "N_dup_names", "T_dup_names", "all_normal_names", "all_tumor_names")
    return(output)

}



