library(dplyr)
library(tibble)

load("./BRCA.exp.rda")
load("./BRCA_dataPrep.rda")

samples_fullNames<-colnames(dataPrep)
samples<-data@colData@listData[["sample"]]
diagnosis<-data@colData@listData[["definition"]]
age_diag<-data@colData@listData[["age_at_index"]]
vital_status<-data@colData@listData[["subtype_vital_status"]]
days_toDeath<-data@colData@listData[["subtype_days_to_death"]]
days_toFollouUp<-data@colData@listData[["days_to_last_follow_up"]]
tumor_stage<-data@colData@listData[["tumor_stage"]]

##outcome=numeric 0 = censuring, 1=dead
##survival=numeric 0 = no survival info, 1=survival info available
##outcome.time=time until end of follow-up, censuring or death in weeks, months or years
metadata<-data.frame(samples_fullNames,diagnosis,age_diag,vital_status,days_toDeath,days_toFollouUp,tumor_stage)
metadata <- metadata %>% add_column(outcome = NA)
metadata <- metadata %>% add_column(survival = NA)

head(metadata)

colnames(metadata)<-c("IDs", "diagnosis", "age", "vital_status", "days_toDeath", "outcome.time", "tumor_stage", "outcome", "survival" )

index=1
for (i in metadata$vital_status) {
    if (is.na(i)) {
        metadata$outcome[index] = NA
        metadata$survival[index] = 0
    }
    else if (i == "Alive"){
        metadata$outcome[index] = 0
        metadata$survival[index] = 1
    }
    else if (i == "Dead"){
        metadata$outcome[index] = 1
        metadata$survival[index] = 1
    }
    index=index+1
}
