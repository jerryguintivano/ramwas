\name{bam.scanBamFile}
\alias{bam.scanBamFile}
\title{Scan a BAM file, filter reads, record read start positions and directions}
\description{
This method scans a BAM file \code{bamfilename}, excludes reads with tag \code{scoretag} below \code{minscore}, 
and records read start positions and directions along with a list of QC measures in an ramwas Rbam object.
}
\usage{
bam.scanBamFile( bamfilename, scoretag = "mapq", minscore = 4 ) 
}
\arguments{
	\item{bamfilename}{
		Name of the BAM file.
	}
	\item{scoretag}{
		The score used to filter the reads. Typically \code{"mapq"} or \code{"AS"}.
	}
	\item{minscore}{
		Alignments with \code{scoretag} below \code{minscore} are excluded.
	}
}
\details{
	.
}
\value{
	The method returns a ramwas \code{Rbam} object with read starts and QC information. 
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

show(rbam$qc);

}
\keyword{ bam }
\keyword{ Rbam }
