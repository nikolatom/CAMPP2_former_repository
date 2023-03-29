## ---- eval = FALSE, echo = TRUE, results='hide'-------------------------------
#  runCampp2(batches=c("tumor_stage","tumor_stage"), prefix="test_CAMPP2",
#  data1=campp2_brca_1, data2=campp2_brca_2, metadata1=campp2_brca_1_meta,
#  metadata2=campp2_brca_2_meta, groups=c("IDs", "diagnosis","IDs", "diagnosis"),
#  technology=c("seq","seq"))

## ---- eval = FALSE------------------------------------------------------------
#  library(CAMPP2)

## ---- eval = FALSE, echo = TRUE, results='hide'-------------------------------
#  load("./data/campp2_brca_1.rda")
#  load("./data/campp2_brca_1_meta.rda")
#  load("./data/campp2_brca_2.rda")
#  load("./data/campp2_brca_2_meta.rda")

## ---- eval = FALSE, echo = TRUE, results='hide'-------------------------------
#  n<-nrow(campp2_brca_1)

## ---- eval = FALSE, echo = TRUE, results='hide'-------------------------------
#  campp2_brca_1_NAs<-apply (campp2_brca_1, 2, function(x) {x[sample( c(1:n), floor(n/10))] <- NA; x} )

## ---- eval = FALSE, echo = TRUE, results='hide'-------------------------------
#  campp2_brca_1_replacedNAs<-ReplaceNAs(data=campp2_brca_1_NAs, pct.NA.row=70, pct.NA.column=80)

## ---- eval = FALSE, echo = TRUE, results='hide'-------------------------------
#  ###Dataset without NA values ("campp2_brca_1_replacedNAs" generated in a previous step) is used as an input.
#  campp2_brca_1_zeroFix<-FixZeros(data=campp2_brca_1_replacedNAs,group=campp2_brca_1_meta$diagnosis, remove.sparse.features=TRUE)

## ---- eval = FALSE, echo = TRUE, results='hide'-------------------------------
#  ###Dataset without zero values ("campp2_brca_1_zeroFix" generated in the previous step) is used as an input.
#  campp2_brca_1_normalized<-NormalizeData(data=campp2_brca_1_zeroFix,group=campp2_brca_1_meta$diagnosis,standardize="TMM",transform="voom",technology="seq")

## ---- eval = FALSE, echo = TRUE, results='hide'-------------------------------
#  ###Normalized data ("campp2_brca_1_normalized" generated in a previous step) are used as an input.
#  campp2_brca_1_batchCorrected<-BatchCorrect(data=campp2_brca_1_normalized,batch=campp2_brca_1_meta$tumor_stage,group=campp2_brca_1_meta$diagnosis,technology="seq")

## ---- eval = FALSE, echo = TRUE, results='hide'-------------------------------
#  ###First 10 features from normalized and batch corrected data "campp2_brca_1_batchCorrected" are used as an input.
#  campp2_brca_1_distributionsFit <- FitDistributions(campp2_brca_1[1:10,])

## ---- eval = FALSE, echo = TRUE, results='hide'-------------------------------
#  ###List of fitted distributions "campp2_brca_1_distributionsFit" and the first
#  ###10 feature counts from example batch corrected data are used as an input.
#  PlotDistributions(campp2_brca_1_batchCorrected[1:10,], campp2_brca_1_distributionsFit)

## ---- eval = FALSE, echo = TRUE, results='hide'-------------------------------
#  ###Normalized and batch corrected data ("campp2_brca_1_batchCorrected" generated in a previous step) are used as an input.
#  meanCounts(campp2_brca_1_batchCorrected, campp2_brca_1_meta$diagnosis)

## ---- eval = FALSE, echo = TRUE, results='hide'-------------------------------
#  ###Normalized and batch corrected data ("campp2_brca_1_batchCorrected" generated in a previous step) are used as an input.
#  PCAPlot(campp2_brca_1_batchCorrected, as.factor(campp2_brca_1_meta$subtype), show.PCA.labels=FALSE, cols=NULL, prefix="test_PCA_plot")

## ---- eval = FALSE, echo = TRUE, results='hide'-------------------------------
#  ###Normalized and batch corrected data ("campp2_brca_1_batchCorrected" generated in a previous step) are used as an input.
#  Kmeans.output<-runKmeans(campp2_brca_1_batchCorrected[1:2000,], num.subsets= NULL, subset.size=NULL, show.PCA.labels = FALSE, colors=NULL, prefix="test", num.km.clusters=NULL, seed=123, pca.scale=FALSE)

## ---- eval = FALSE, echo = TRUE, results='hide'-------------------------------
#  ###Sub-sampled (the first 2000 genes) normalized and batch corrected data ("campp2_brca_1_batchCorrected" generated before) are used as an input.
#  EstimateKmeans(t(campp2_brca_1_batchCorrected[1:2000,]))

## ---- eval = FALSE, echo = TRUE, results='hide'-------------------------------
#  campp2_brca_1_DEA<-RunDEA(data=campp2_brca_1_normalized,
#  metadata=campp2_brca_1_meta,
#  group=campp2_brca_1_meta$subtype, prefix="test",
#  block=campp2_brca_1_meta$subtype, batch=campp2_brca_1_meta$age,
#  covarDEA = c("tumor_stage"), cutoff.logFC=1, cutoff.FDR=0.01)

## ---- eval = FALSE, echo = TRUE, results='hide'-------------------------------
#  DEA_all_comparisons <- DEAFeatureApply(data = campp2_brca_1_normalized,
#  design.matrix = campp2_brca_1_DEA$DEA.design.matrix, contrast.matrix =
#  campp2_brca_1_DEA$DEA.contrast.matrix, cutoff.logFC =1, cutoff.FDR =0.01,
#  block = campp2_brca_1_meta$subtype, vector = FALSE)

## ---- eval = FALSE, echo = TRUE, results='hide'-------------------------------
#  DEA_one_comparison <-
#  DEAFeature(contrast.matrix = campp2_brca_1_DEA$DEA.contrast.matrix[,1],
#  data = campp2_brca_1_normalized,
#  design.matrix = campp2_brca_1_DEA$DEA.design.matrix,
#  cutoff.logFC =1, cutoff.FDR =0.01, block = campp2_brca_1_meta$subtype)

## ---- eval = FALSE, echo = TRUE, results='hide'-------------------------------
#  campp2_brca_1_DEA_out<-ExportDEA(res.DEA = campp2_brca_1_DEA$res.DEA,
#  prefix="test")

## ---- eval = FALSE, echo = TRUE, results='hide'-------------------------------
#  RunDEAVisuals(campp2_brca_1_DEA_HUGO, cutoff.FDR = 0.01, cutoff.logFC = 1,
#  n.labeled.features = 15, control.group= "healthy", prefix="test_DEA_visuals")

## ---- eval = FALSE, echo = TRUE, results='hide'-------------------------------
#  MakeVolcano(campp2_brca_1_DEA_HUGO, prefix = "test_volcano",
#  cutoff.logFC = 1, cutoff.FDR = 0.01, n.labeled.features = 15)
#  

## ---- eval = FALSE, echo = TRUE, results='hide'-------------------------------
#  MakeUpset("testUpSet",campp2_brca_1_DEA_HUGO_features_per_group,
#  label="test_UpSet",label.vertical.position=16.5,
#  y.axis.by=300,set.size.by=250)

## ---- eval = FALSE, echo = TRUE, results='hide'-------------------------------
#  RunHeatmap(campp2_brca_1_batchCorrected, campp2_brca_1_DEA_HUGO,
#  campp2_brca_1_meta$subtype, heatmap.size=30, viridis.palette="turbo",
#  plot.heatmap="DEA", "test_heatmap")

## ---- eval = FALSE, echo = TRUE, results='hide'-------------------------------
#  MakeHeatmap(data=campp2_brca_1_batchCorrected, group=campp2_brca_1_meta$subtype,
#  prefix="test_heatmap", data.type = "expression/abundance")
#  

## ---- eval = FALSE, echo = TRUE, results='hide'-------------------------------
#  ###Normalized and batch corrected data are used as an input for regression.
#  ##run regression using small batch corrected dataset - no double cross validation
#  runLASSO(data=campp2_brca_1_batchCorrected, group=campp2_brca_1_meta$diagnosis, alpha=0.5, min.coef=0, nfolds=10, prefix="test")
#  ##run regression using large batch corrected dataset suitable for double cross validation
#  ###create a large dataset
#  campp2_test_data_LASSO<-cbind(campp2_brca_1,campp2_brca_2)
#  campp2_test_metadata_LASSO<-rbind(campp2_brca_1_meta, campp2_brca_2_meta)
#  campp2_test_data_LASSO_replaceNAs<-ReplaceNAs(data=campp2_test_data_LASSO)
#  campp2_test_data_LASSO_zeroFix<-FixZeros(data=campp2_test_data_LASSO_replaceNAs,group=campp2_test_metadata_LASSO$diagnosis, remove.sparse.features=TRUE)
#  campp2_test_data_LASSO_normalized<-NormalizeData(data=campp2_test_data_LASSO_zeroFix,group=campp2_test_metadata_LASSO$diagnosis,standardize="TMM",transform="voom",technology="seq")
#  campp2_test_data_LASSO_batchCorrected<-BatchCorrect(data=campp2_test_data_LASSO_normalized,batch=campp2_test_metadata_LASSO$tumor_stage,group=campp2_test_metadata_LASSO$diagnosis,technology="seq")
#  ###run lasso on a large dataset
#  runLASSO(data=campp2_test_data_LASSO_batchCorrected, group=campp2_test_metadata_LASSO$diagnosis, alpha=0.5, min.coef=0, nfolds=10, prefix="test")

## ---- eval = FALSE, echo = TRUE, results='hide'-------------------------------
#  ###Normalized and batch corrected data are used as an input for regression.
#  LASSOFeature(seed=123, data=campp2_brca_1_batchCorrected,
#  group=campp2_brca_1_meta$diagnosis, alpha=0.5, validation=TRUE,
#  min.coef = 0, nfolds=10)

## ---- eval = FALSE, echo = TRUE, results='hide'-------------------------------
#  ##here we prepare DEA results compatible (group="diagnosis) with
#  EN/LASS/Ridge regression results:
#  campp2_brca_1_DEA_diagnosis<-RunDEA(data=campp2_brca_1_normalized,
#  metadata=campp2_brca_1_meta,
#  group=campp2_brca_1_meta$diagnosis, prefix="test", batch=campp2_brca_1_meta$age,
#  covarDEA = c("tumor_stage"), cutoff.logFC=1, cutoff.FDR=0.01)
#  ##here we run RunDEA_LASSO_consensus function:
#  campp2_brca_1_DEA_LASSO_consensus <- RunDEA_LASSO_consensus(
#  campp2_brca_1_DEA_diagnosis$DEA.out,
#  campp2_brca_1_LASSO, group=campp2_brca_1_meta$diagnosis,
#  viridis.palette="turbo", "test")

## ---- eval = FALSE, echo = TRUE, results='hide'-------------------------------
#  campp2_brca_1_run_rf <-
#  RunRF(data = campp2_brca_1_batchCorrected,
#  group = campp2_brca_1_meta$diagnosis,
#  split.size = 5,
#  test.train.ratio = 0.25,
#  num.trees.init = 5000,
#  num.trees.iterat = 2000)

## ---- eval = FALSE, echo = TRUE, results='hide'-------------------------------
#  campp2_brca_1_rf_apply <- RFApply(data = campp2_brca_1_batchCorrected,
#  group = campp2_brca_1_meta$diagnosis,
#  validation = TRUE,
#  test.train.ratio = 0.25,
#  num.trees.init = 5000,
#  num.trees.iterat = 2000)

## ---- eval = FALSE, echo = TRUE, results='hide'-------------------------------
#  campp2_brca_1_forest_features <- ForestFeatures(seed = 123,
#  data = campp2_brca_1_batchCorrected,
#  group = campp2_brca_1_meta$diagnosis,
#  validation = TRUE,
#  test.train.ratio = 0.25,
#  num.trees.init = 5000,
#  num.trees.iterat = 2000)

