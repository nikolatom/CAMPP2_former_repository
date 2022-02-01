calc_stat<- function(gene_counts,sample){
    gene_counts<-as.data.frame(gene_counts)
    replicate<-gene_counts %>% select(contains(sample))
    sumRep<-""
    for (i in 1:ncol(replicate)){
        sumRep[i]<-sum(replicate[,i])
    }
    sumRep<-as.integer(sumRep)
    sumALL<-sum(sumRep)

    norm<-NULL
    normRep<-NULL
    for (i in 1:ncol(replicate)){
        norm[i]<-sumALL/sumRep[i]
        normRep[i]<-data.frame(norm[i]*(replicate[,i]))
    }
    summary(as.data.frame(normRep))
}

calc_stat(dataPrep,"A0DB-01")
calc_stat(dataPrep,"A0DC-01")
calc_stat(dataPrep,"A13D-01")
calc_stat(dataPrep,"A13E-01")
calc_stat(dataPrep,"A26E-01")
calc_stat(dataPrep,"A26J-01")

"A0DB-01" 30, 3.9, 12 (7.7 : 1 : 3.1)
"A0DC-01" 8, 28 (1 : 3.5)
"A13D-01" 2, 12, 90 (1 : 6 : 45)
"A13E-01" 30.9, 2.3, 11 (13.4 : 1 : 4.8)
"A26E-01" 38, 9, 3.7 (10.2 : 2.4 : 1)
"A26J-01" 3.8, 40, 6 (1: 10.5 : 1.6)

