
#Is there a better way how to create global variables?
parseArguments2 <- function(data, sdata, metadata, smetadata, groups, technology, batches, datacheck, standardize, transform, plotmds, plotheatmap, kmeans, signif, colors, prefix, corr, lasso, WGCNA, cutoffWGCNA, survival, covar, stratify, survplot, PPint, GenemiRint){
    print(technology)
    print(groups)
    my_list<-list("data"=data, "sdata"=sdata, "groups"=groups, "technology"=technology, "prefix"=prefix)
    return(my_list)
#    print(sdata)
#   print(prefix)
}
