## ----loadPackages, echo=FALSE, warning=FALSE, message=FALSE--------------
library(knitr)
library(pander)
suppressPackageStartupMessages(library(ramwas))
panderOptions("digits", 3)
# opts_chunk$set(eval=FALSE)
# dr = "D:/temp/"

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
writeLines(con = paste0(dr,"/CpG_chromosome_names.txt"), text = chrnames)

## ----setseed3, echo=FALSE------------------------------------------------
set.seed(18090212)

## ----fillDataMat---------------------------------------------------------
fmm = fm.create(paste0(dr,"/Coverage"), nrow = nsamples, ncol = nvariables)
fms = fm.create(paste0(dr,"/SNPs"), nrow = nsamples, ncol = nvariables,
                size = 1, type = "integer")

# Row names of the matrices are set to sample names
rownames(fmm) = as.character(covariates$sample)
rownames(fms) = as.character(covariates$sample)

# The matrices are filled, 2000 variables at a time
byrows = 2000
for( i in seq_len(nvariables/byrows) ){ # i=1
    ind = (1:byrows) + byrows*(i-1)

    snps = rbinom(n = byrows * nsamples, size = 2, prob = 0.2)
    dim(snps) = c(nsamples, byrows)
    fms[,ind] = snps

    slice = double(nsamples*byrows)
    dim(slice) = c(nsamples, byrows)
    slice[,  1:225] = slice[,  1:225] + covariates$sex / 50 / sd(covariates$sex)
    slice[,101:116] = slice[,101:116] + covariates$age / 16 / sd(covariates$age)
    slice = slice +
            ((as.integer(factor(covariates$batch))+i) %% 3) / 200 +
            snps/1.5 +
            runif(nsamples*byrows)/2
    fmm[,ind] = slice;
}
close(fms)
close(fmm)

## ----paramMWAS, warning=FALSE, message=FALSE-----------------------------
param = ramwasParameters(
    dircoveragenorm = dr,
    covariates = covariates,
    modelcovariates = "batch",
    modeloutcome = "sex",
    toppvthreshold = 20,
    fileSNPs = "SNPs"
)

## ----threads, echo=FALSE-------------------------------------------------
# Bioconductor requires limit of 2 parallel jobs
param$cputhreads = 2

## ----MWAS, message=FALSE-------------------------------------------------
ramwas5MWAS(param)

## ----plotQQ1, echo=FALSE, warning=FALSE, message=FALSE, eval=TRUE, fig.width=6, fig.height=6----
pfull = parameterPreprocess(param)
fm = fm.open( paste0(pfull$dirmwas, "/Stats_and_pvalues") )
pv = fm[,3]
close(fm)
qqPlotFast(pv)
title(pfull$qqplottitle)

## ----SNPs, message=FALSE-------------------------------------------------
ramwasSNPs(param)

## ----plotQQ2, echo=FALSE, warning=FALSE, message=FALSE, eval=TRUE, fig.width=6, fig.height=6----
pfull = parameterPreprocess(param)
fm = fm.open( paste0(pfull$dirSNPs, "/Stats_and_pvalues") )
pv = fm[,3]
close(fm)
qqPlotFast(pv)
title("QQ-plot for MWAS with correction for SNPs")

## ----topPvMWAS-----------------------------------------------------------
# Get the directory with testing results
pfull = parameterPreprocess(param)
toptbl = read.table(
                paste0(pfull$dirmwas,"/Top_tests.txt"),
                header = TRUE,
                sep = "\t")
pander(head(toptbl,10))
toptbl = read.table(
                paste0(pfull$dirSNPs,"/Top_tests.txt"),
                header = TRUE,
                sep = "\t")
pander(head(toptbl,10))

