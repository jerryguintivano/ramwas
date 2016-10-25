## ----loadPackages, echo=FALSE, warning=FALSE, message=FALSE--------------
suppressPackageStartupMessages(library(ramwas))
suppressPackageStartupMessages(library(BSgenome.Ecoli.NCBI.20080805))

## ----cpgsFromGenome, warning=FALSE, message=FALSE------------------------
library(ramwas)
library(BSgenome.Ecoli.NCBI.20080805)
cpgset = getCpGsetCG(BSgenome.Ecoli.NCBI.20080805)
# First 10 CpGs in NC_008253:
print(cpgset$NC_008253[1:10])

## ----getCpGsetALL1, eval=FALSE, warning=FALSE, message=FALSE-------------
#  library(BSgenome.Hsapiens.UCSC.hg19)
#  library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
#  genome = injectSNPs(Hsapiens, "SNPlocs.Hsapiens.dbSNP144.GRCh37")
#  cpgset = getCpGsetALL(genome)
#  # Number of CpGs with all SNPs injected in autosomes
#  sum(sapply(cpgset[1:22], length))

## ----echo=FALSE----------------------------------------------------------
42841152

## ----getCpGsetALL2, eval=FALSE-------------------------------------------
#  # Do for all chromosomes
#  genome[['chr22']] =
#      injectSNPsMAF(
#          gensequence = BSGenome[['chr22']],
#          frqcount = 'count_ALL_chr22.txt',
#          MAF = 0.01)
#  
#  # Find the CpGs
#  cpgset = getCpGsetALL(genome)

## ----save1, eval=FALSE---------------------------------------------------
#  saveRDS(file = 'My_cpgset.rds', object = cpgset)

## ----insilicoFASTQ, eval=FALSE-------------------------------------------
#  # Do for all chromosomes
#  insilicoFASTQ(
#      con="chr1.fastq.gz",
#      gensequence = BSGenome[['chr1']],
#      fraglength=75)

## ----RaMWAS, eval=FALSE--------------------------------------------------
#  library(ramwas)
#  chrset = paste0('chr',1:22)
#  targetcov = 75
#  covtolerance = 10
#  
#  param = list(
#  	dirproject = '.',
#  	dirbam = './bams',
#  	dirfilter = TRUE,
#  	bamnames = chrset,
#  	bam2sample = list(all_samples = chrset),
#  	scoretag = "AS",
#  	minscore = 100,
#  	minfragmentsize = targetcov,
#  	maxfragmentsize = targetcov,
#  	minavgcpgcoverage = 0,
#  	minnonzerosamples = 0,
#  	# filecpgset - file with the CpG set being QC-ed
#  	filecpgset = filecpgset
#  )
#  param1 = parameterPreprocess(param)
#  ramwas1scanBams(param)
#  ramwas3NormalizedCoverage(param)

## ----filter, eval=FALSE--------------------------------------------------
#  # Preprocess parameters to learn the location of coverage matrix
#  param1 = parameterPreprocess(param)
#  
#  # Load the coverage matrix (vector)
#  cover = fm.load( paste0(param1$dircoveragenorm, '/Coverage'))
#  
#  # split the coverage by chromosomes
#  # `cpgset` - the CpG set being QC-ed
#  fac = rep(seq_along(cpgset), times = sapply(cpgset, length))
#  levels(fac) = names(cpgset)
#  class(fac) = 'factor'
#  cover = split(cover, fac)
#  	
#  # filter CpGs on each chromosome by the coverage
#  cpgsetQC = cpgset
#  for( i in seq_along(cpgset) ){
#      keep = (cover[[i]] >= (targetcov - covtolerance)) &
#             (cover[[i]] <= (targetcov + covtolerance))
#      cpgsetQC[[i]] = cpgset[[i]][ keep ]
#  }

## ----save2, eval=FALSE---------------------------------------------------
#  saveRDS(file = 'My_cpgset_QC.rds', object = cpgsetQC)

