## ----global_options, include=FALSE---------------------------------------
#getwd()
#knitr::opts_chunk$set(fig.align="center", fig.retina=1)
knitr::opts_chunk$set(fig.retina=1)
library(ramwas)

## ----loadCgGset----------------------------------------------------------
filename = system.file("extdata", "bigQC.rds", package = "ramwas");
qc = readRDS(filename)$qc;
# show(names(qc))


## ----nbams---------------------------------------------------------------
cat("N BAMs:", qc$nbams)

## ----reads.total---------------------------------------------------------
cat("Reads total:", qc$reads.total)

## ----reads.aligned-------------------------------------------------------
{
 cat("Reads aligned:", qc$reads.aligned, "\n")
 cat("This is ", qc$reads.aligned / qc$reads.total * 100, 
     "% of all reads", sep="")
}

## ----reads.recorded------------------------------------------------------
{
 cat("Reads recorded:",qc$reads.recorded,"\n")
 cat("This is ", qc$reads.recorded / qc$reads.aligned * 100, 
     "% of aligned reads", sep="")
}

## ----frwrev--------------------------------------------------------------
{
 cat("Reads on forward strand:", qc$frwrev[1],"\n")
 cat("Reads on reverse strand:", qc$frwrev[2],"\n")
 cat("Fraction of reads on forward strand:", qcmean(qc$frwrev), "\n")
}

## ----hist.score1, fig.width=8--------------------------------------------
{
 cat("Average alignment score:", qcmean(qc$hist.score1), "\n")
 cat("Average alignment score, no filter:", qcmean(qc$bf.hist.score1), "\n")
 par(mfrow=c(1,2)) 
 plot(qc$hist.score1)
 plot(qc$bf.hist.score1)
}

## ----hist.length.matched, fig.width=8------------------------------------
{
 cat("Average aligned length:", qcmean(qc$hist.length.matched), "\n")
 cat("Average aligned length, no filter:", 
     qcmean(qc$bf.hist.length.matched), "\n")
 par(mfrow = c(1,2)) 
 plot(qc$hist.length.matched)
 plot(qc$bf.hist.length.matched)
}

## ----hist.edit.dist1-----------------------------------------------------
{
 cat("Average edit distance:", qcmean(qc$hist.edit.dist1), "\n")
 cat("Average edit distance, no filter:", qcmean(qc$bf.hist.edit.dist1), "\n")
 par(mfrow = c(1,2)) 
 plot(qc$hist.edit.dist1)
 plot(qc$bf.hist.edit.dist1)
}

## ----reads.recorded.no.repeats-------------------------------------------
{
 cat("Reads without duplicates:", qc$reads.recorded.no.repeats, "\n")
 cat("This is ", qc$reads.recorded.no.repeats / qc$reads.recorded * 100, 
     "% of aligned reads", "\n", sep="")
 cat("Fraction of reads on forward strand (with    duplicates):", 
     qcmean(qc$frwrev), "\n")
 cat("Fraction of reads on forward strand (without duplicates):", 
     qcmean(qc$frwrev.no.repeats), "\n")
}

## ----cnt.nonCpG.reads----------------------------------------------------
{
 cat("Non-CpG reads:", qc$cnt.nonCpG.reads[1], "\n")
 cat("This is ", qcmean(qc$cnt.nonCpG.reads)*100, "% of recorded reads", sep="")
}

## ----avg.cpg.coverage----------------------------------------------------
{
 cat("Summed across", qc$nbams, "bams", "\n")
 cat("Average     CpG coverage:", qc$avg.cpg.coverage,   "\n")
 cat("Average non-CpG coverage:", qc$avg.noncpg.coverage,"\n")
 cat("Enrichment ratio:", qc$avg.cpg.coverage / qc$avg.noncpg.coverage)
}

## ----hist.isolated.dist1-------------------------------------------------
plot(qc$hist.isolated.dist1)

## ----chrXY---------------------------------------------------------------
{
 cat("ChrX reads: ", qc$chrX.count[1], ", which is ", 
     qcmean(qc$chrX.count)*100, "% of total", sep="", "\n")
 cat("ChrX reads: ", qc$chrY.count[1], ", which is ", 
     qcmean(qc$chrY.count)*100, "% of total", sep="", "\n")
}

## ----avg.coverage.by.density---------------------------------------------
{
	cat("Highest coverage is observed at CpG density of", 
	    qcmean(qc$avg.coverage.by.density)^2);
	plot(qc$avg.coverage.by.density)
}

