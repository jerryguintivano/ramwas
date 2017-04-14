## ----loadKnitr, echo=FALSE-----------------------------------------------
# library("knitr")
# opts_chunk$set(eval=FALSE)
library(pander)
panderOptions("digits", 3)

## ----eval=FALSE----------------------------------------------------------
#  param = ramwasParameters(
#      dirproject = ".",
#      dirbam = "bams",
#      filebamlist = "bam_list.txt",
#      filecpgset = "Simulated_chromosome.rds",
#      cputhreads = 2,
#      scoretag = "MAPQ",
#      minscore = 4,
#      minfragmentsize = 50,
#      maxfragmentsize = 250,
#      filecovariates = "covariates.txt",
#      modelcovariates = NULL,
#      modeloutcome = "age",
#      modelPCs = 0,
#      toppvthreshold = 1e-5,
#      cvnfolds = 10,
#      mmalpha = 0,
#      mmncpgs = c(5,10,50,100,500,1000,5000,10000)
#  )

## ----eval=FALSE----------------------------------------------------------
#  ### R parameter file
#  dirbam = "/ramwas_project/bams/"
#  dirproject = "/ramwas_project/"
#  filebamlist = "/ramwas_project/000_list_of_files.txt"
#  scoretag = "AS"
#  minscore = 100
#  
#  ### platform dependent part
#  if(.Platform$OS.type == "windows"){
#      filecpgset="C:/RaMWAS/CpG_set/cpgset_hg19_SNPS_at_MAF_0.05.rds"
#  } else {
#      filecpgset="/computing_cluster/ramwas/cpgset_hg19_SNPS_at_MAF_0.05.rds"
#  }

## ----eval=FALSE----------------------------------------------------------
#  batch1/b1sample1.bam
#  batch1/b1sample2.bam
#  batch2/b2sample1.bam
#  batch2/b2sample2.bam
#  batch2/b2sample3.bam
#  batch4/sample4.bam

## ----eval=FALSE----------------------------------------------------------
#  batch1/sample1.bam
#  batch1/sample2.bam
#  batch2/sample1.bam
#  batch2/sample2.bam

## ----bam2sample----------------------------------------------------------
bam2sample = list(
    sample1 = c("bam1","bam2","bam3"),
    sample2 = "sample2"
)

## ----CpGsetExample-------------------------------------------------------
cpgset = list( chr1 = c(12L, 57L, 123L),
               chr2 = c(45L, 95L, 99L, 111L),
               chr3 = c(22L, 40L, 199L, 211L) )

## ----marts, eval=FALSE---------------------------------------------------
#  library(biomaRt)
#  library(ramwas)
#  
#  # First pick a host.
#  bihost = "grch37.ensembl.org"
#  
#  # First we list databases
#  listOfMarts = listMarts(host = bihost)
#  pander(head(listOfMarts, 10))
#  
#  # Pick a database
#  bimart = "ENSEMBL_MART_ENSEMBL"
#  
#  # Connect to the database
#  mart = useMart(biomart = bimart, host = bihost)
#  
#  # List the data sets in the database
#  listOfDatasets = listDatasets(mart = mart)
#  pander(head(listOfDatasets, 10))
#  
#  # Pisk a data set
#  bidataset = "hsapiens_gene_ensembl"
#  
#  # Connect to the data set
#  mart = useMart(biomart = bimart, dataset = bidataset, host = bihost)
#  
#  # List the attributes
#  listOfAttributes = listAttributes(mart)
#  pander(head(listOfAttributes, 10))
#  
#  # Pick attributes
#  biattributes = c("hgnc_symbol", "entrezgene", "strand")
#  
#  listOfFilters = listFilters(mart)
#  pander(head(listOfFilters, 20))
#  
#  # Pick a filter
#  bifilters = list(with_hgnc_trans_name=TRUE)
#  
#  # Test a location
#  chr = "chr1"
#  pos =  15975530
#  param = ramwasParameters(
#      bihost = bihost,
#      bimart = bimart,
#      bidataset = bidataset,
#      biattributes = biattributes,
#      bifilters = bifilters,
#      biflank = 0);
#  
#  anno = ramwasAnnotateLocations(param, chr, pos);
#  pander(anno)

