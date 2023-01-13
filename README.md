Cancer Systems Biology, Section of Bioinformatics, Department of Health and Technology, Technical University of Denmark, 2800, Lyngby, Copenhagen, Denmark

Cancer Structural Biology Group, Danish Cancer Society Research Center, Strandboulevarden 49, 2100, Copenhagen, Denmark

# CAncer-bioMarker-Prediction-Pipeline - CAMPP2  #
<br/>

### Installation instructions

#### Create and activate conda environment with R and devtools in your project directory

If you don't have access to conda please see the Miniconda installer page (https://docs.conda.io/en/latest/miniconda.html) on instructions on how to install Miniconda.

Remember to have installed BiocManager upfront before proceeding with the installation of CAMPP2.   
install.packages("BiocManager")


Once you have installed it:

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
The data for testing the functions and workflow includes 2 BRCA datasets (campp2_brca_1, campp2_brca_2) and the associated metadata (campp2_brca_1_meta, campp2_brca_2_meta). Each dataset is represented by raw read counts (10000 genes) for 30 samples: 20 tumours which are divided into 4 subtypes (each subtype has 5 samples), and 10 normals. Metadata includes information about diagnosis, age, vital status, days to death, outcome time, tumor stage, subtype, outcome and survival.

Both, raw read counts and metadata were extracted from TCGA BReast CAncer dataset (TCGA-BRCA) Level 3 data. 

Test data (gene counts and metadata for 2 data sets) are integrated into CAMPP2 package and accessible as campp2_brca_1, campp2_brca_2, campp2_brca_1_meta, campp2_brca_2_meta variables once the package is installed. Results from intermediate steps (not described here) are also integrated and used as an input for running the examples of the functions. R data object files (.rda) are available on https://github.com/ELELAB/CAMPP2/tree/main/data. 


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
`runCampp2(batches=c("tumor_stage","tumor_stage"),prefix="test_CAMPP2", data1=campp2_brca_1, data2=campp2_brca_2, metadata1=campp2_brca_1_meta,metadata2=campp2_brca_2_meta, groups=c("IDs", "diagnosis","IDs", "diagnosis"), technology=c("seq","seq"))`
<br/>

For testing the functions, you can consider the code present in `campp2_example_Run.R` (git repository). <br/>

For more details, see help page of the main function `runCampp2`.
