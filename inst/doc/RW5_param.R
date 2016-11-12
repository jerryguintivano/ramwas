## ----eval=FALSE----------------------------------------------------------
#  param = list(
#      dirbam = "/ramwas_project/bams/",
#      dirproject = "/ramwas_project/",
#      filebamlist = "/ramwas_project/000_list_of_files.txt",
#      scoretag = "AS",
#      minscore = 100,
#      cputhreads = 8,
#      filecpgset = "/ramwas/cpgsets/cpgset_hg19_SNPS_at_MAF_0.05.rds",
#      maxrepeats = 3,
#      maxfragmentsize=200,
#      minfragmentsize=50,
#      bamnames = NULL,
#      filebam2sample = "/ramwas_project/bam2sample.txt",
#      minavgcpgcoverage = 0.3,
#      minnonzerosamples = 0.3,
#      covfile = "cov.txt",
#      chrkeep = 1:22
#  );

## ----eval=FALSE----------------------------------------------------------
#  ### R parameter file
#  dirbam = "/ramwas_project/bams/"
#  dirproject = "/ramwas_project/"
#  filebamlist = "/ramwas_project/000_list_of_files.txt"
#  scoretag = "AS"
#  minscore = 100
#  
#  ### platform dependent part
#  if(.Platform$OS.type == "windows") {
#      filecpgset='C:/RaMWAS/CpG_set/cpgset_hg19_SNPS_at_MAF_0.05.rds'
#  } else {
#      filecpgset='/computing_cluster/ramwas/cpgset_hg19_SNPS_at_MAF_0.05.rds'
#  }
#  
#  

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

## ----CpGsetExample-------------------------------------------------------
cpgset = list( chr1 = c(12L, 57L, 123L),
               chr2 = c(45L, 95L, 99L, 111L),
               chr3 = c(22L, 40L, 199L, 211L) );

## ----eval=FALSE----------------------------------------------------------
#  sample1=b1sample1,b2sample1
#  sample2=b1sample2,b2sample2
#  sample3=b2sample3
#  sample4

