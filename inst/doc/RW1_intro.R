## ----loadKnitr, echo=FALSE-----------------------------------------------
library("knitr")
# opts_chunk$set(eval=FALSE)
library(pander)
panderOptions("digits", 3)

## ----install, eval=FALSE-------------------------------------------------
#  ## try http:// if https:// URLs are not supported
#  source("http://bioconductor.org/biocLite.R")
#  biocLite("ramwas")

## ----loadIt, eval=FALSE--------------------------------------------------
#  library(ramwas) # Loads the package
#  browseVignettes("ramwas") # Opens vignettes
#  help(package="ramwas") # Lists package functions

## ----loadPackages, echo=FALSE, warning=FALSE, message=FALSE--------------
suppressPackageStartupMessages(library(ramwas))
# dr = "D:/temp"

## ----generateData, warning=FALSE-----------------------------------------
library(ramwas)
dr = paste0(tempdir(), "/simulated_project")
ramwas0createArtificialData(dir = dr,
                            verbose = FALSE,
                            nreads = 100e3,
                            ncpgs = 25e3)
# This is the project directory
cat(dr)

## ----parameters----------------------------------------------------------
param = ramwasParameters(
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
    modelcovariates = "sex",
    modeloutcome = "age",
    modelPCs = 0,
    toppvthreshold = 1e-5,
    bihost = "grch37.ensembl.org",
    bimart = "ENSEMBL_MART_ENSEMBL",
    bidataset = "hsapiens_gene_ensembl",
    biattributes = c("hgnc_symbol","entrezgene","strand"),
    bifilters = list(with_hgnc_trans_name=TRUE),
    biflank = 0,
    cvnfolds = 5,
    mmalpha = 0,
    mmncpgs = c(5,10,50,100,500,1000,2000,3000)
)

## ----scan-bams, warning=FALSE, message=FALSE-----------------------------
ramwas1scanBams(param)

## ----plotACbD, echo=FALSE, warning=FALSE, message=FALSE, fig.width=6, fig.height=6----
pfull = parameterPreprocess(param)
qc = readRDS(paste0(pfull$dirrqc, "/BAM007.qc.rds"))
plot(qc$qc$avg.coverage.by.density)

## ----collectQC1, warning=FALSE, message=FALSE----------------------------
ramwas2collectqc(param)

## ----plotFSD, echo=FALSE, warning=FALSE, message=FALSE, fig.width=6, fig.height=6----
qc = readRDS(paste0(pfull$dirqc, "/summary_total/qclist.rds"))
frdata = qc$total$hist.isolated.dist1
estimate = as.double(readLines(
    paste0(pfull$dirproject,"/Fragment_size_distribution.txt")))
ramwas:::plotFragmentSizeDistributionEstimate(frdata, estimate)

## ----normCoverage99, warning=FALSE, message=FALSE------------------------
ramwas3normalizedCoverage(param)

## ----pca99, warning=FALSE, message=FALSE---------------------------------
ramwas4PCA(param)

## ----plotPCA, echo=FALSE, warning=FALSE, message=FALSE, fig.width=6, fig.height=6----
e = readRDS(paste0(pfull$dirpca, "/eigen.rds"))
ramwas:::plotPCvalues(e$values)
ramwas:::plotPCvectors(e,2)
rm(e)

## ----tablePCAcr, echo=FALSE, warning=FALSE, message=FALSE----------------
tblcr = read.table(paste0(pfull$dirpca, "/PC_vs_covs_corr.txt"),
                 header = TRUE,
                 sep = "\t")
pander(head(tblcr, 3))

## ----tablePCApv, echo=FALSE, warning=FALSE, message=FALSE----------------
tblpv = read.table(paste0(pfull$dirpca, "/PC_vs_covs_pvalue.txt"),
                 header = TRUE,
                 sep = "\t")
pander(head(tblpv, 3))

## ----mwas99, warning=FALSE, message=FALSE--------------------------------
ramwas5MWAS(param)

## ----tableMWAS, echo=FALSE, warning=FALSE, message=FALSE, fig.width=6, fig.height=6----
fm = fm.open( paste0(pfull$dirmwas, "/Stats_and_pvalues") )
pv = fm[,3]
close(fm)
qqPlotFast(pv)
title(pfull$qqplottitle)

## ----anno, warning=FALSE, message=FALSE, eval=FALSE----------------------
#  ramwas6annotateTopFindings(param)

## ----CV, warning=FALSE, message=FALSE------------------------------------
ramwas7riskScoreCV(param)

## ----plotCV1, echo=FALSE, warning=FALSE, message=FALSE, fig.width=6, fig.height=6----
cv = readRDS( paste0(pfull$dircv, "/rds/CpGs=000050_alpha=0.000000.rds") )
ramwas:::plotPrediction(param = pfull,
                        outcome = cv$outcome,
                        forecast = cv$forecast,
                        cpgs2use = 50,
                        main = "Prediction success (EN on coverage)")

## ----plotCV2, echo=FALSE, warning=FALSE, message=FALSE, fig.width=6, fig.height=6----
ramwas7CplotByNCpGs(param)
cl = readRDS(sprintf("%s/rds/cor_data_alpha=%f.rds",
                    pfull$dircv,
                    pfull$mmalpha))
ramwas:::plotCVcors(cl, pfull)

## ----topPvMWAS-----------------------------------------------------------
# Get the directory with testing results
pfull = parameterPreprocess(param)
toptbl = read.table(
                paste0(pfull$dirmwas,"/Top_tests.txt"),
                header = TRUE,
                sep = "\t")
pander(head(toptbl, 5))

## ----getLocation---------------------------------------------------------
chr = toptbl$chr[1]
position = toptbl$position[1]
cat("Top Finding is at:", chr, "-", position, "\n")

## ----getdata-------------------------------------------------------------
datavec = getDataByLocation(param, chr, position)
testres = getTestsByLocation(param, chr, position)
pander(testres)

## ----lm------------------------------------------------------------------
outcome = pfull$covariates[[pfull$modeloutcome]]
cvrt = pfull$covariates[[pfull$modelcovariates]]
variable = datavec$matrix
model = lm( outcome ~ variable + cvrt)
pander(summary(model)$coefficients)

## ----dirlocations--------------------------------------------------------
pfull = parameterPreprocess(param)
# Here lies coverage matrix
pfull$dircoveragenorm
# Here are PCA files
pfull$dirpca
# Here are MWAS files
pfull$dirmwas
# Here are multi-marker cross-validation files
pfull$dircv

## ----clean---------------------------------------------------------------
unlink(paste0(dr,"/*"), recursive=TRUE)

## ----version-------------------------------------------------------------
sessionInfo()

