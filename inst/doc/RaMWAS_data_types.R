## ----global_options, include=FALSE---------------------------------------
#knitr::opts_chunk$set(fig.align='center', fig.retina=1)

## ----loadCgGset----------------------------------------------------------
filename = system.file("extdata", "cpgset.rds", package = "ramwas");
cpgset = readRDS(filename);
show(cpgset)


## ----CpGsetExample-------------------------------------------------------
cpgset = list( chr1 = c(12L, 57L, 123L),
               chr2 = c(45L, 95L, 99L, 111L),
               chr3 = c(22L, 40L, 199L, 211L) );

## ----scanBam-------------------------------------------------------------
suppressPackageStartupMessages(library(ramwas))
filename = system.file("extdata", "small.bam", package = "ramwas");
rbam = bam.scanBamFile(filename, scoretag = "AS", minscore = 80 ) 

## ----rbamReads-----------------------------------------------------------
head(names(rbam$startsfwd))
rbam$startsfwd$chr2
rbam$startsrev$chr1


## ----loadSampleQC--------------------------------------------------------
filename = system.file("extdata", "qc_sample.rds", package = "ramwas");
rbam = readRDS(filename);

## ----reads.total---------------------------------------------------------
rbam$qc$reads.total

## ----reads.aligned-------------------------------------------------------
rbam$qc$reads.aligned

## ----reads.recorded------------------------------------------------------
rbam$qc$reads.recorded
# sum(sapply(rbam$startsfwd,length)) + sum(sapply(rbam$startsrev,length))

## ----frwrev--------------------------------------------------------------
rbam$qc$frwrev
# c(sum(sapply(rbam$startsfwd,length)), sum(sapply(rbam$startsrev,length)))

## ----hist.length.matched, fig.align='center'-----------------------------
barplot(rbam$qc$hist.length.matched/1e6, width = 1, space = 0, col = "blue", main = "Distribution of the length of the aligned part of the reads", xaxs="i", yaxs="i", ylab = "count, millions");
at = seq(0,length(rbam$qc$hist.length.matched)+10,10);
at[1] = 1;
axis(1,at = at-0.5, labels = at)

## ----hist.edit.dist1-----------------------------------------------------
barplot(rbam$qc$hist.edit.dist1/1e6, width = 1, space = 0, col = "blue", main = "Distribution of edit distance", xaxs="i", yaxs="i", ylab = "count, millions");
at = seq(0,length(rbam$qc$hist.edit.dist1)+10,1);
axis(1,at = at+0.5, labels = at)

## ----hist.score1---------------------------------------------------------
barplot(rbam$qc$hist.score1/1e6, width = 1, space = 0, col = "blue", main = "Distribution of read scores", xaxs="i", yaxs="i", ylab = "count, millions");
at = seq(0,length(rbam$qc$hist.score1)+10,20);
axis(1,at = at+0.5, labels = at)

