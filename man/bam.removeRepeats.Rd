\name{bam.removeRepeats}
\alias{bam.removeRepeats}
\title{Remove highly repeated reads from an Rbam object.}
\description{
As a side effect of DNA amplification some reads can be sequenced multiple times as if they were independent reads.
These excessive reads do not represent real signal and should be filtered out.

To filter out such reads the function \code{bam.removeRepeats} finds reads that start at the same position and have the same direction (forward or reverse strand) and removes the reads is excess of \code{maxrep} (maximum repeat) quantity.
}
\usage{
bam.removeRepeats(rbam, maxrep)
}
\arguments{
	\item{rbam}{
		Rbam ramwas object.
	}
	\item{maxrep}{
		Maximum number of reads allowed to start at the same position and have the same direction.
	}
}
\details{
	.
}
\value{
	The method returns a ramwas \code{Rbam} object with excess repeat reads removed.
	The object has two extra QC measures added, \code{frwrev.no.repeats} and \code{reads.recorded.no.repeats}.
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

rbam2 = bam.removeRepeats(rbam, 3)

cat('After removing repeats.','\n');
cat('Reads on forward strand:',rbam2$qc$frwrev.no.repeats[1],'\n');
cat('Reads on reverse strand:',rbam2$qc$frwrev.no.repeats[2],'\n');

}
\keyword{ bam }
\keyword{ Rbam }
