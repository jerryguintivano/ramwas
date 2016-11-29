## ----switcher, echo=FALSE------------------------------------------------
library("knitr")
# opts_chunk$set(eval=FALSE)

## ----install, eval=FALSE-------------------------------------------------
#  ## try http:// if https:// URLs are not supported
#  source("http://bioconductor.org/biocLite.R")
#  biocLite("ramwas")

## ----loadIt, eval=FALSE--------------------------------------------------
#  library(ramwas) # Loads the package
#  browseVignettes("ramwas") # Opens vignettes
#  help(package="ramwas") # Lists package functions

## ----loadPackages, echo=FALSE, warning=FALSE, message=FALSE, eval=TRUE----
suppressPackageStartupMessages(library(ramwas))

## ----generateData, warning=FALSE-----------------------------------------
library(ramwas)
dr = paste0(tempdir(), "/simulated_project")
ramwas0createArtificialData(dir = dr, verbose = FALSE)
# This is the project directory
cat(dr)

## ----parameters, eval=TRUE-----------------------------------------------
param = list(
    dirproject = dr,
    dirbam = "bams",
    filebamlist = "bam_list.txt",
    filecpgset = "Simulated_chromosome.rds",
    cputhreads = 2,
    scoretag = "MAPQ",
    minscore = 4,
    minfragmentsize = 50,
    maxfragmentsize = 250,
    minavgcpgcoverage = 0.3,
    minnonzerosamples = 0.3,
    filecovariates = "covariates.txt",
    modelcovariates = NULL,
    modeloutcome = "age",
    modelPCs = 0,
    toppvthreshold = 1e-5,
    bihost = "grch37.ensembl.org",
    bimart = "ENSEMBL_MART_ENSEMBL",
    bidataset = "hsapiens_gene_ensembl",
    biattributes = c("hgnc_symbol","entrezgene","strand"),
    bifilters = list(with_hgnc_transcript_name=TRUE),
    biflank = 0,
    cvnfolds = 10,
    mmalpha = 0,
    mmncpgs = c(5,10,50,100,500,1000,5000,10000)
)

## ----scan-bams, warning=FALSE, message=FALSE-----------------------------
ramwas1scanBams(param)

## ----collectQC1, warning=FALSE, message=FALSE----------------------------
ramwas2collectqc(param)

## ----normCoverage99, warning=FALSE, message=FALSE------------------------
ramwas3normalizedCoverage(param)

## ----pca99, warning=FALSE, message=FALSE---------------------------------
ramwas4PCA(param)

## ----mwas99, warning=FALSE, message=FALSE--------------------------------
ramwas5MWAS(param)

## ----anno, warning=FALSE, message=FALSE----------------------------------
ramwas6annotateTopFindings(param)

## ----parameters7, warning=FALSE, message=FALSE---------------------------
ramwas7riskScoreCV(param)

## ----dirlocations, eval=TRUE---------------------------------------------
fullparam = parameterPreprocess(param)
# Here lies coverage matrix
fullparam$dircoveragenorm
# Here are PCA files
fullparam$dirpca
# Here are MWAS files
fullparam$dirmwas
# Here are multi-marker cross-validation files
fullparam$dircv

## ----clean---------------------------------------------------------------
unlink(paste0(dr,"/*"), recursive=TRUE)

## ----version, eval=TRUE--------------------------------------------------
sessionInfo()

