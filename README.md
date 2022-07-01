Cancer Systems Biology, Section of Bioinformatics, Department of Health and Technology, Technical University of Denmark, 2800, Lyngby, Copenhagen, Denmark

Cancer Structural Biology Group, Danish Cancer Society Research Center, Strandboulevarden 49, 2100, Copenhagen, Denmark

# CAncer-bioMarker-Prediction-Pipeline - CAMPP2  #
<br/>

### Installation instructions

#### Create and activate conda environment with R and devtools in your project directory
`conda create --prefix env -c conda-forge r-base=4.1.3 r-devtools` <br/>
`conda activate ./env`


#### Install CAMPP2 from Github
open R console: <br/>
`R` <br/>


load devtools library <br/>
`library(devtools)`

install CAMPP2 from private Github repository (using your personal token and commit ID)
`devtools::install_github(repo = "ELELAB/CAMPP2",auth_token="<your personal token to github>")`

### Example data
The data to test the functions and workflow includes: 2 BRCA datasets (campp2_brca_1, campp2_brca_2) and the associated metadata (campp2_brca_1_meta, campp2_brca_2_meta). Each dataset is represented by raw read counts for 25 samples (15 tumours, 10 normals in the 1st dataset; 18 tumours and 7 normals in the 2nd dataset). Metadata includes information about diagnosis, age, vital status, days to death, outcome time, tumor stage, outcome and survival.

Both, raw read counts and metadata were extracted from TCGA BReast CAncer dataset (TCGA-BRCA) Level 3 data. 

Test data are integrated into CAMPP2 package and accessible as campp2_brca_1, campp2_brca_2, campp2_brca_1_meta, campp2_brca_2_meta variables once the package is installed. R data object files (.rda) are also available on https://github.com/ELELAB/CAMPP2/tree/main/data. 


### Example run
Default settings can be executed in R using:

`library(CAMPP2)`

Test data are already part of the CAMPP2 package so user doesn't need to download them. In case you want to load .rda objects manually from the cloned repository, you can use this code:

`load("./data/campp2_brca_1.rda")` 
<br/>
`load("./data/campp2_brca_1_meta.rda")` 
<br/>
`load("./data/campp2_brca_2.rda")` 
<br/>
`load("./data/campp2_brca_2_meta.rda")` 
<br/>
<br/>
Default workflow could be run using this command: 
<br/>
`runCampp2(batches=c("tumor_stage","tumor_stage"),prefix="test_CAMPP2", data1=dataset1, data2=dataset2, metadata1=metadata1,metadata2=metadata2, groups=c("IDs", "diagnosis","IDs", "diagnosis"), technology=c("seq","seq"))`
<br/>

For testing the functions, you can consider the code present in `campp2_example_Run.R` (git repository). <br/>

For more details, see help page of the main function `runCampp2`.
