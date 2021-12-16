library(stringr)
library(dplyr) 
library(tibble)

setwd("")
featurecounts<-read.table(file = "all.featureCounts.tsv", sep = "\t",)
table1<-subset (featurecounts, select = -c(2:6))

header<-str_extract(table1[1,], "\\d+[_]+\\d+[.T.T02.S2]+") 
header[1]<-"Geneid"

table1[1,]<-header

names(table1) <- table1[1,]
table1 <- table1[-1,]

samples_to_remove<-c("11021_000000.T.T02.S2",
                     "21191_000000.T.T02.S2",
                     "23944_000000.T.T02.S2",
                     "19603_000001.T.T02.S2",
                     "22901_000000.T.T02.S2",
                     "22243_000000.T.T02.S2",
                     "22901_000001.T.T02.S2",
                     "24409_000000.T.T02.S2",
                     "24261_000000.T.T02.S2",
                     "21057_000000.T.T02.S2",
                     "21094_000001.T.T02.S2",
                     "21371_000000.T.T02.S2",
                     "22025_000001.T.T02.S2",
                     "22197_000000.T.T02.S2",
                     "23250_000001.T.T02.S2",
                     "19603_000000.T.T02.S2",
                     "21135_000000.T.T02.S2")

table_final<-select(table1,-starts_with(samples_to_remove))

write.table(table_final, file = "campp_input_featureCounts_gencode.txt",row.names = FALSE, sep = "\t")

