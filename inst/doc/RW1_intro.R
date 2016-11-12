## ----switcher, echo=FALSE------------------------------------------------
library('knitr')
opts_chunk$set(eval=FALSE)

## ----install, eval=FALSE-------------------------------------------------
#  ## try http:// if https:// URLs are not supported
#  source("http://bioconductor.org/biocLite.R")
#  biocLite("ramwas")

## ----loadIt, eval=FALSE--------------------------------------------------
#  library(ramwas) # Loads the package
#  browseVignettes('ramwas') # Opens vignettes
#  help(package="ramwas") # Lists package functions

## ----loadPackages, echo=FALSE, warning=FALSE, message=FALSE--------------
#  suppressPackageStartupMessages(library(ramwas))

## ----generateData, warning=FALSE-----------------------------------------
#  library(ramwas)
#  dr = paste0(tempdir(), "/simulated_project");
#  
#  ramwas0createArtificialData(dir = dr,
#                              nsamples = 20,
#                              nreads = 1e6,
#                              ncpgs = 250e3,
#                              verbose=FALSE)
#  # This is the project directory
#  cat(dr)

## ----parameters0---------------------------------------------------------
#  param = list(
#      dirproject = dr,
#      dirbam = "bams",
#      filebamlist = "bam_list.txt",
#      filecpgset = "Simulated_chromosome.rds"
#  )

## ----cputhreads, echo=FALSE----------------------------------------------
#  param$cputhreads = 2

## ----parameters1---------------------------------------------------------
#  param$scoretag = "AS"
#  param$minscore = 4

## ----parameters2---------------------------------------------------------
#  param$minfragmentsize = 50
#  param$maxfragmentsize = 250

## ----scan-bams, warning=FALSE, message=FALSE-----------------------------
#  ramwas1scanBams(param)

## ----parameters3---------------------------------------------------------
#  param$filebam2sample = "bam_list.txt"

## ----collectQC1, warning=FALSE, message=FALSE----------------------------
#  ramwas2collectqc(param)

## ----parameters4---------------------------------------------------------
#  param$minavgcpgcoverage = 0.3
#  param$minnonzerosamples = 0.3

## ----normCoverage99, warning=FALSE, message=FALSE------------------------
#  ramwas3NormalizedCoverage(param)

## ----parameters5---------------------------------------------------------
#  param$filecovariates = "covariates.txt"
#  param$modelcovariates = NULL

## ----pca99, warning=FALSE, message=FALSE---------------------------------
#  ramwas4PCA(param)

## ----parameters6---------------------------------------------------------
#  param$modeloutcome = "casecontrol"
#  param$modelPCs = 0

## ----mwas99, warning=FALSE, message=FALSE--------------------------------
#  ramwas5MWAS(param)

## ----mwas88, warning=FALSE, message=FALSE--------------------------------
#  param$modeloutcome = "age"
#  param$modelPCs = 0
#  ramwas5MWAS(param)

## ----parameters7, warning=FALSE, message=FALSE---------------------------
#  param$modeloutcome = "age"
#  param$modelPCs = 0
#  param$cvnfolds = 10;
#  ramwas6crossValidation(param)

## ----cv1, warning=FALSE, message=FALSE-----------------------------------
#  param$mmalpha = 0;
#  param$mmncpgs = 50;
#  ramwas7multiMarker(param)

## ----cv2, warning=FALSE, message=FALSE-----------------------------------
#  param$mmalpha = 0;
#  param$mmncpgs = 10
#  ramwas7multiMarker(param)

## ----clean---------------------------------------------------------------
#  unlink(dr, recursive=TRUE)

## ----version, eval=TRUE--------------------------------------------------
sessionInfo()

