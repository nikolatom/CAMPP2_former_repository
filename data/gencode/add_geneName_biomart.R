require(biomaRt)

ensembl104 <- useEnsembl(biomart = 'genes', 
                        dataset = 'hsapiens_gene_ensembl',
                        version = 104)

#listAttributes(ensembl104)
setwd("/home/projects/dtu_00011/data/shared_projects/campp2_ALL/campp/CAMPP_analysis/data/CAMPP_input/gencode")
data <- read.table("campp_input_featureCounts_gencode.txt", header = FALSE)
data[,1] <- sub('\\.[0-9]*$', '', data[,1])
data[1,1] <- "Gene stable ID"
colnames(data) <- data[1,]
data<-data[-1,]

ensembl_id <- as.data.frame(data[,1])
colnames(ensembl_id) <- "Gene stable ID"

gene_names <- getBM(
  attributes = c('ensembl_gene_id', 'external_gene_name'),
  mart = ensembl104,
  values = ensembl_id,
  filters = 'ensembl_gene_id',
  uniqueRows = TRUE,
  bmHeader = T)

data_merged <- merge(gene_names, data, by="Gene stable ID", all = TRUE)
data_merged
write.table(data_merged,file="campp_input_featureCounts_gencode_hgnc.txt")
