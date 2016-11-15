## ----eval=FALSE----------------------------------------------------------
#  param = list(
#      dirproject = ".",
#      dirbam = "bams",
#      filebamlist = "bam_list.txt",
#      filecpgset = "Simulated_chromosome.rds",
#      cputhreads = 2,
#      scoretag = "MAPQ",
#      minscore = 4,
#      minfragmentsize = 50,
#      maxfragmentsize = 250,
#      minavgcpgcoverage = 0.3,
#      minnonzerosamples = 0.3,
#      filecovariates = "covariates.txt",
#      modelcovariates = NULL,
#      modeloutcome = "age",
#      modelPCs = 0,
#      toppvthreshold = 1e-5,
#      bidataset = "hsapiens_gene_ensembl",
#      biattributes = c("hgnc_symbol","entrezgene","strand"),
#      bifilters = list(with_hgnc_transcript_name=TRUE),
#      biflank = 0,
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

