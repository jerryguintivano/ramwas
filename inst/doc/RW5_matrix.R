## ----loadPackages, echo=FALSE, warning=FALSE, message=FALSE--------------
library(knitr)
library(pander)
suppressPackageStartupMessages(library(ramwas))
panderOptions("digits", 3)
# opts_chunk$set(eval=FALSE)
# dr = "D:/temp/"
nsamples = 200

## ----generateData--------------------------------------------------------
library(ramwas)

# work in a temporary directory
dr = paste0(tempdir(), "/simulated_matrix_data")
dir.create(dr, showWarnings = FALSE)
cat(dr,"\n")

## ----dims----------------------------------------------------------------
nsamples = 200
nvariables = 100000

## ----setseed1, echo=FALSE------------------------------------------------
set.seed(18090212)

## ----genCovar------------------------------------------------------------
covariates = data.frame(
    sample = paste0("Sample_",seq_len(nsamples)),
    sex = seq_len(nsamples) %% 2,
    age = runif(nsamples, min = 20, max = 80),
    batch = paste0("batch",(seq_len(nsamples) %% 3))
)
pander(head(covariates))

## ----setseed2, echo=FALSE------------------------------------------------
set.seed(18090212)

## ----genLocs-------------------------------------------------------------
temp = cumsum(sample(20e7 / nvariables, nvariables, replace = TRUE) + 0)
chr      = as.integer(temp %/% 1e7) + 1L
position = as.integer(temp %% 1e7)

locmat = cbind(chr = chr, position = position)
chrnames = paste0("chr", 1:22)
pander(head(locmat))

## ----locSave-------------------------------------------------------------
fmloc = fm.create.from.matrix(
            filenamebase = paste0(dr,"/CpG_locations"),
            mat = locmat)
close(fmloc)
writeLines(con = paste0(dr,"/CpG_chromosome_names.txt"),
           text = chrnames)

## ----setseed3, echo=FALSE------------------------------------------------
set.seed(18090212)

## ----fillDataMat---------------------------------------------------------
fm = fm.create(paste0(dr,"/Coverage"), nrow = nsamples, ncol = nvariables)

# Row names of the matrix are set to sample names
rownames(fm) = as.character(covariates$sample)

# The matrix is filled, 2000 variables at a time
byrows = 2000
for( i in seq_len(nvariables/byrows) ){ # i=1
    slice = matrix(runif(nsamples*byrows), nrow = nsamples, ncol = byrows)
    slice[,  1:225] = slice[,  1:225] + covariates$sex / 30 / sd(covariates$sex)
    slice[,101:116] = slice[,101:116] + covariates$age / 10 / sd(covariates$age)
    slice = slice + ((as.integer(factor(covariates$batch))+i) %% 3) / 40
    fm[,(1:byrows) + byrows*(i-1)] = slice
}
close(fm)

## ----param1--------------------------------------------------------------
param = ramwasParameters(
    dircoveragenorm = dr,
    covariates = covariates,
    modelcovariates = NULL
)

## ----threads, echo=FALSE-------------------------------------------------
# Bioconductor requires limit of 2 parallel jobs
param$cputhreads = 2

## ----pcaNULL, warning=FALSE, message=FALSE-------------------------------
ramwas4PCA(param)

## ----plotPCA, echo=FALSE, warning=FALSE, message=FALSE, fig.width=6, fig.height=6----
pfull = parameterPreprocess(param)
e = readRDS(paste0(pfull$dirpca, "/eigen.rds"))
ramwas:::plotPCvalues(e$values)
ramwas:::plotPCvectors(e,1)
ramwas:::plotPCvectors(e,2)
ramwas:::plotPCvectors(e,3)
ramwas:::plotPCvectors(e,4)
rm(e)

## ----topCorNULL----------------------------------------------------------
# Get the directory with PCA results
pfull = parameterPreprocess(param)
tblcr = read.table(paste0(pfull$dirpca, "/PC_vs_covs_corr.txt"),
                 header = TRUE,
                 sep = "\t")
pander(head(tblcr, 5))

## ----topPvNULL-----------------------------------------------------------
pfull = parameterPreprocess(param)
tblpv = read.table(paste0(pfull$dirpca, "/PC_vs_covs_pvalue.txt"),
                 header = TRUE,
                 sep = "\t")
pander(head(tblpv, 5))

## ----pcaBatch, warning=FALSE, message=FALSE------------------------------
param$modelcovariates = "batch"

ramwas4PCA(param)

## ----topPvBatch----------------------------------------------------------
# Get the directory with PCA results
pfull = parameterPreprocess(param)
tblpv = read.table(paste0(pfull$dirpca, "/PC_vs_covs_pvalue.txt"),
                 header = TRUE,
                 sep = "\t")
pander(head(tblpv, 5))

## ----paramGWAS, warning=FALSE, message=FALSE-----------------------------
param$modelcovariates = "batch"
param$modeloutcome = "sex"
param$toppvthreshold = 20

ramwas5MWAS(param)

## ----tableMWAS, echo=FALSE, warning=FALSE, message=FALSE, fig.width=6, fig.height=6----
pfull = parameterPreprocess(param)
fm = fm.open( paste0(pfull$dirmwas, "/Stats_and_pvalues") )
pv = fm[,3]
close(fm)
qqPlotFast(pv)
title(pfull$qqplottitle)

## ----topPvMWAS-----------------------------------------------------------
# Get the directory with testing results
pfull = parameterPreprocess(param)
toptbl = read.table(
                paste0(pfull$dirmwas,"/Top_tests.txt"),
                header = TRUE,
                sep = "\t")
pander(head(toptbl, 5))

## ----clean---------------------------------------------------------------
unlink(paste0(dr,"/*"), recursive=TRUE)

## ----version-------------------------------------------------------------
sessionInfo()

