#' @title RegExp extraction of gene counts
#' @description This function extracts gene counts based on sample name using regular expressions. This might be useful e.g. when working with TCGA datasets.
#' @param gene_count_matrix gene count matrix to be processed
#' @param normal_sampleRegExp regular expression to extract normal samples (based on the sample name); default value fits TCGA data (https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/)
#' @param tumor_sampleRegExp regular expression to extract tumor samples
#' @export
#' @import dplyr
#' @seealso
#' @return As an output, list of 1) read counts from paired normal samples ($normal_paired), 2) read counts from paired tumor samples ($tumor_paired), 3) counts from all normal samples ($normal), 4)  counts from all tumor samples ($tumor), 5) sample names of duplicated normal samples ($N_dup_names), 6) sample names of duplicated tumor samples ($T_dup_names), 7) sample names from all normal samples ($all_normal_names), 8) sample names from all tumor samples ($all_tumor_names) is generated.
#' @examples \dontrun{
#'     load("~/opt/campp_bioconducor/CAMPP2/data/BRCA/BRCA_dataPrep.rda")
#'     normal_tumor_counts<-extract_tumor_normal_samples(dataPrep)
#'     normal_counts_paired<-as.data.frame(normal_tumor_counts$normal_paired)
#'     tumor_counts_paired<-as.data.frame(normal_tumor_counts$tumor_paired) #including technical replicates
#'     normal_counts<-as.data.frame(normal_tumor_counts$normal)
#'     tumor_counts<-as.data.frame(normal_tumor_counts$tumor)
#'     normal_duplicates<-normal_tumor_counts$N_dup_names
#'     tumor_duplicates<-normal_tumor_counts$T_dup_names
#'     all_normal_names<-normal_tumor_counts$all_normal_names
#'     all_tumor_names<-normal_tumor_counts$all_tumor_names
#'     write.table(tumor_counts_paired, file="tumor_counts_paired.txt", sep = "\t")
#'     write.table(normal_counts_paired, file="normal_counts_paired.txt", sep = "\t")
#'     write.table(tumor_counts, file="tumor_counts.txt", sep = "\t")
#'     write.table(normal_counts, file="normal_counts.txt", sep = "\t")
#'     write.table(tumor_duplicates, file="tumor_duplicates.txt", sep = "\t")
#'     write.table(normal_duplicates, file="normal_duplicates.txt", sep = "\t")
#' }

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



