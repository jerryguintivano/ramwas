\name{bam.QCs}
\alias{bam.count.nonCpG.reads}
\alias{bam.hist.isolated.distances}
\title{Functions for quality control measures for Rbam files.}
\description{
The functions listed here calculate quality control measures in addition to those calculated by \code{\link{bam.scanBamFile}} and \code{\link{bam.removeRepeats}}.

The function \code{bam.count.nonCpG.reads} counts the number of reads that would not cover any CpGs if they were \code{distance} bp in length. With \code{distance} set large enough, the fraction of reads not covering any CpGs should be low (<2\%). 
The function adds \code{$cnt.nonCpG.reads} QC measure which is a vector of two numbers.
First is the number of reads not covering CpGs. The second is the total number of reads for convenient calculation of the ratio.
The reads on chromosomes not present in the CpG set are not included in calculations.

The function \code{bam.hist.isolated.distances} takes reads covering isolated CpGs at distance \code{distance} and calculated the distribution of the distances from the read start positions to the corresponding isolated CpGs.
The function adds \code{$hist.isolated.dist1} QC measure which contains the distribution information. The first element contains the number of reads starting at the isolated CpG. The element 2 contains the number of reads starting at distance 1 from the CpG (C on the corresponding strand).

}
\usage{
bam.count.nonCpG.reads(rbam, cpgset, distance)
bam.hist.isolated.distances(rbam, isocpgset, distance)
}
\arguments{
	\item{rbam}{
		Rbam ramwas object.
	}
	\item{cpgset}{
		Set of CpG locations.
	}
	\item{isocpgset}{
		Set of locations of isolated CpGs.
	}	
	\item{distance}{
		Minimum distance to any CpGs in \code{cpgset}.
	}
}
\details{
	.
}
\value{
	The functions return a ramwas \code{Rbam} object with the corresponding QC added.
	
	
}
\references{
	.
}
\seealso{
	.
}
\note{
	.
}
\author{
	Andrey Shabalin \email{ashabalin@vcu.edu}
}
\examples{

filename = system.file("extdata", "small.bam", package = "ramwas");

rbam = bam.scanBamFile( filename, scoretag = "AS", minscore = 80 ) 

data(toycpgset);
show(toycpgset)

rbam2 = bam.count.nonCpG.reads(rbam, toycpgset, 50)
cat( rbam2$qc$cnt.nonCpG.reads[1], 'of',
     rbam2$qc$cnt.nonCpG.reads[2], 'reads are not covering CpGs', '\n');

rbam3 = bam.hist.isolated.distances(rbam, toycpgset, 50)

}
\keyword{ bam }
\keyword{ Rbam }
