\name{QCs}
\alias{qcmean}
\alias{qcmean.NULL}
\alias{qcmean.qcChrX}
\alias{qcmean.qcChrY}
\alias{qcmean.qcCoverageByDensity}
\alias{qcmean.qcEditDist}
\alias{qcmean.qcEditDistBF}
\alias{qcmean.qcFrwrev}
\alias{qcmean.qcHistScore}
\alias{qcmean.qcHistScoreBF}
\alias{qcmean.qcIsoDist}
\alias{qcmean.qcLengthMatched}
\alias{qcmean.qcLengthMatchedBF}
\alias{qcmean.qcNonCpGreads}
\alias{plot.qcCoverageByDensity}
\alias{plot.qcEditDist}
\alias{plot.qcEditDistBF}
\alias{plot.qcHistScore}
\alias{plot.qcHistScoreBF}
\alias{plot.qcIsoDist}
\alias{plot.qcLengthMatched}
\alias{plot.qcLengthMatchedBF}
\title{
Quality Control Measures
}
\description{
RaMWAS calculates a number of QC measures for each
BAM and saves them in R .rds files.

For full description of the QC measures and the ploting options
run\cr \code{vignette("RW3_BAM_QCs")}
}
\usage{
qcmean(x)
\method{qcmean}{NULL}(x)
\method{qcmean}{qcChrX}(x)
\method{qcmean}{qcChrY}(x)
\method{qcmean}{qcCoverageByDensity}(x)
\method{qcmean}{qcEditDist}(x)
\method{qcmean}{qcEditDistBF}(x)
\method{qcmean}{qcFrwrev}(x)
\method{qcmean}{qcHistScore}(x)
\method{qcmean}{qcHistScoreBF}(x)
\method{qcmean}{qcIsoDist}(x)
\method{qcmean}{qcLengthMatched}(x)
\method{qcmean}{qcLengthMatchedBF}(x)
\method{qcmean}{qcNonCpGreads}(x)

\method{plot}{qcHistScore}(x, samplename="", xstep = 25, ...)
\method{plot}{qcHistScoreBF}(x, samplename="", xstep = 25, ...)
\method{plot}{qcEditDist}(x, samplename="", xstep = 5, ...)
\method{plot}{qcEditDistBF}(x, samplename="", xstep = 5, ...)
\method{plot}{qcLengthMatched}(x, samplename="", xstep = 25, ...)
\method{plot}{qcLengthMatchedBF}(x, samplename="", xstep = 25, ...)
\method{plot}{qcIsoDist}(x, samplename="", xstep = 25, ...)
\method{plot}{qcCoverageByDensity}(x, samplename="", ...)
}
\arguments{
    \item{x}{
        The QC object. See the examples below.
    }
    \item{samplename}{
        Name of the sample for plot title.
    }
    \item{xstep}{
        The distance between x axis ticks.
    }
    \item{...}{
        Parameters passed to the underlying
        \code{\link[graphics]{plot}} or
        \code{\link[graphics]{barplot}} function.
    }
}

\value{
    Function \code{qcmean} returns one value summary of most QC measures.
    Run \code{vignette("RW3_BAM_QCs")} for description of values returned by it.

}
\author{
    Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\seealso{
    See vignettes: \code{browseVignettes("ramwas")}.
}
\examples{
# Load QC data from a sample project
filename = system.file("extdata", "bigQC.rds", package = "ramwas")
qc = readRDS(filename)$qc


## The number of BAM files
cat("N BAMs:", qc$nbams)

## Total number of reads in the BAM file(s)
cat("Reads total:", qc$reads.total)

## Number of reads aligned to the reference genome
cat("Reads aligned:", qc$reads.aligned, "\\n")
cat("This is ", qc$reads.aligned / qc$reads.total * 100,
    "\% of all reads", sep="")

## Number of reads that passed minimum score filter and are recorded
cat("Reads recorded:",qc$reads.recorded,"\\n")
cat("This is ", qc$reads.recorded / qc$reads.aligned * 100,
    "\% of aligned reads", sep="")

## Number of recorded reads aligned to
## the forward and reverse strands respectively
cat("Reads on forward strand:", qc$frwrev[1],"\\n")
cat("Reads on reverse strand:", qc$frwrev[2],"\\n")
cat("Fraction of reads on forward strand:", qcmean(qc$frwrev), "\\n")

## Distribution of the read scores
cat("Average alignment score:", qcmean(qc$hist.score1), "\\n")
cat("Average alignment score, no filter:", qcmean(qc$bf.hist.score1), "\\n")
par(mfrow=c(1,2))
plot(qc$hist.score1)
plot(qc$bf.hist.score1)

## Distribution of the length of the aligned part of the reads
cat("Average aligned length:", qcmean(qc$hist.length.matched), "\\n")
cat("Average aligned length, no filter:",
    qcmean(qc$bf.hist.length.matched), "\\n")
par(mfrow = c(1,2))
plot(qc$hist.length.matched)
plot(qc$bf.hist.length.matched)

## Distribution of edit distance between
## the aligned part of the read and the reference genome
cat("Average edit distance:", qcmean(qc$hist.edit.dist1), "\\n")
cat("Average edit distance, no filter:", qcmean(qc$bf.hist.edit.dist1), "\\n")
par(mfrow = c(1,2))
plot(qc$hist.edit.dist1)
plot(qc$bf.hist.edit.dist1)

## Number of reads after removal of duplicate reads
cat("Reads without duplicates:", qc$reads.recorded.no.repeats, "\\n")
cat("This is ", qc$reads.recorded.no.repeats / qc$reads.recorded * 100,
    "\% of aligned reads", "\\n", sep="")
cat("Fraction of reads on forward strand (with    duplicates):",
    qcmean(qc$frwrev), "\\n")
cat("Fraction of reads on forward strand (without duplicates):",
    qcmean(qc$frwrev.no.repeats), "\\n")

## Number of reads away from CpGs
cat("Non-CpG reads:", qc$cnt.nonCpG.reads[1], "\\n")
cat("This is ", qcmean(qc$cnt.nonCpG.reads)*100, "\% of recorded reads",
    sep="")

## Average coverage of CpGs and non-CpGs
cat("Summed across", qc$nbams, "bams", "\\n")
cat("Average     CpG coverage:", qc$avg.cpg.coverage,   "\\n")
cat("Average non-CpG coverage:", qc$avg.noncpg.coverage,"\\n")
cat("Enrichment ratio:", qc$avg.cpg.coverage / qc$avg.noncpg.coverage)

## Coverage around isolated CpGs
plot(qc$hist.isolated.dist1)

## Fraction of reads from chrX and chrY
cat("ChrX reads: ", qc$chrX.count[1], ", which is ",
    qcmean(qc$chrX.count)*100, "\% of total", sep="", "\\n")
cat("ChrX reads: ", qc$chrY.count[1], ", which is ",
    qcmean(qc$chrY.count)*100, "\% of total", sep="", "\\n")

## Coverage vs. CpG density
cat("Highest coverage is observed at CpG density of",
    qcmean(qc$avg.coverage.by.density)^2)
plot(qc$avg.coverage.by.density)
}
