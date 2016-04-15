\name{estimateFragmentSizeDistribution}
\alias{estimateFragmentSizeDistribution}
\title{Estimate distribution of fragment sizes based on distribution of reads around isolated CpG.}
\description{
The function fits a logistic curve to the distribution of distances from read starts to isolated CpGs.
}
\usage{
estimateFragmentSizeDistribution(hist.isolated.distances, seqLength)
}
\arguments{
	\item{hist.isolated.distances}{
		Distribution of distances from read starts to isolated CpGs.
		Usually calculated with \code{\link{bam.hist.isolated.distances}} and accumulated accross all samples in the project.		}
	\item{seqLength}{
		The minimum size of DNA fragments (in bp). Usually set to the sequenced length of the fragments.
	}
}
\details{
	If length of \code{hist.isolated.distances} is less or equal to \code{seqLength}
}
\value{
	Returns a vector of probabilities that a fragment is at least that long or longer.
	Naturally, first \code{seqLength} elements are equal to 1. The remaining monotonically decay to zero.
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
0
}
\keyword{ bam }
\keyword{ Rbam }
