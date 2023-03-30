Cancer Systems Biology, Section of Bioinformatics, Department of Health and Technology, Technical University of Denmark, 2800, Lyngby, Copenhagen, Denmark

Cancer Structural Biology Group, Danish Cancer Society Research Center, Strandboulevarden 49, 2100, Copenhagen, Denmark

# CAncer-bioMarker-Prediction-Pipeline - CAMPP2  #
<br/>

### Installation instructions

#### Installation from BioConductor 

To install the CAMPP2 package from Bioconductor, you need to follow these steps:

1. Install the Bioconductor package manager, BiocManager, if you don't have it already. 
You can do this by running the following command in your R console:


`if (!requireNamespace("BiocManager", quietly = TRUE))` <br/>
    `install.packages("BiocManager")` 


2. Install the CAMPP2 package using BiocManager by running the following command:

`BiocManager::install("CAMPP2")` 


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

### CAMPP2 development environment

The CAMPP2 development environment is made available as Docker image, defined
in a Dockerfile in this repository. It uses the Bioconductor devel Docker image
as a base, and further downloads and installs the requirements for CAMPP2.
In order to use it, please install [Docker](https://www.docker.com) if you
don't have it already, and follow these instructions:

1. Clone the CAMPP2 GitHub repository in a local folder:

```
git clone https://www.github.com/ELELAB/CAMPP2.git
cd CAMPP2
```

2. Build the container image:

```
docker build --pull -t campp2:devel-20230206 ./docker
```

We recommend tagging the container with a date consistent with the day it's been
created, since it might contain different package versions depending
on the date, as done in this example.

Notice that this needs to be done only once, as well as every time you intend
to upgrade your image to account for the latest Bioconductor devel docker
image (see below)

3. Run a container from the image you just created:

In order to access the development environment via RStudio on browser, run:

```
docker run \
    --rm \
    -v /path/to/my/CAMPP2:/home/rstudio/CAMPP2 \
    -p 8787:8787 \
    -e PASSWORD=campp2 \
    campp2:devel-20230206
```

Open your web browser and head to http://localhost:8787. Log in using `rstudio`
as username and `campp2` as password. You should be able to access the CAMPP2
development folder from the `Rstudio` interface.

If you prefer running the R command line prompt (e.g. using `R`, it works in
the same exacy way with `Rscript`):

```
docker run \
    --rm \
    -it \
    -v /path/to/my/CAMPP2:/home/rstudio/CAMPP2 \
    campp2:devel-20230206 \
    R
```

if you want to access the `bash` command line prompt of the container you
can instead run:

```
docker run \
    --rm \
    -it \
    -v /path/to/my/CAMPP2:/home/rstudio/CAMPP2 \
    campp2:devel-20230206 \
    bash
```

Notice that in these commands:
  - the `/path/to/my/CAMPP2` should be replaced by the actual absolute path to the
  CAMPP2 development folder, depending on where it is located in your computer
  - the container name should be changed to the actual name of the container
  that you built in step 2. If you're unsure what the container name or tag is, run
  `docker image list`, you should be able to find both (a container is specified as
  `NAME:TAG`)

We recommend to rebuild your image every few weeks to account for changes in
the original BioConductor development image and in the CRAN/Bioconductor
repositories. To do so, just repeat step 2 and adjust step 3 according to
your new tag.

#### Running on an Apple Silicon Macs

If you are running the container on Apple Silicon (i.e. M1 and M2 at the time
of writing) we currently recommend running the container under x86 emulation
using Rosetta2. In order to do so,

  1. Install Rosetta2 by running on your terminal:

```
/usr/sbin/softwareupdate --install-rosetta
```

  2. Turn on the "Use Rosetta for x86/amd64 emulation on Apple Silicon" option
in Docker, currently located in Settings, Features under development

Further, we recommend adding the `--platform linux/amd64` option to the
`docker run` command lines above, for example:

```
docker run \
    --rm \
    -it \
    -v /path/to/my/CAMPP2:/home/rstudio/CAMPP2 \
    --platform linux/amd64 \
    campp2:devel-20230206 \
    bash
```

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
