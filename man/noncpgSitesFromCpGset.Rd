\name{noncpgSitesFromCpGset}
\alias{noncpgSitesFromCpGset}
\alias{isocpgSitesFromCpGset}
\title{Find a set of genomic locations away from CpGs / a set of CpGs away from other CpGs.}
\description{
Several QC measures requires a set of locations away from any CpGs.
Several others use a set of CpGs located away from any other CpGs.

The function \code{noncpgSitesFromCpGset} creates a set of locations away from all CpGs from a set of CpG locations.
The set of produced by picking a middle location between each 2 consecutive CpGs at least \code{2 * distance} apart.

The function \code{isocpgSitesFromCpGset} creates a set of CpGs away from all other CpGs from a set of CpG locations.
}
\usage{
noncpgSitesFromCpGset(cpgset, distance)
isocpgSitesFromCpGset(cpgset, distance)
}
\arguments{
	\item{cpgset}{
		Set of CpG locations.
	}
	\item{distance}{
		Minimum distance to any CpGs / other CpGs in \code{cpgset}.
	}
}
\details{
	.
}
\value{
	Both functions return a set of locations for the same set of chromosomes as in \code{cpgset}.

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

data(toycpgset)

show(toycpgset)

toynoncpgset = noncpgSitesFromCpGset(toycpgset, 10)

show(toynoncpgset)

toyisocpgset = isocpgSitesFromCpGset(toycpgset, 10)

show(toyisocpgset)

### Test minimum distance for a non-CpG to fit in.

cpgset = lapply(1:10, function(x){return(c(1,1+x,1+2*x))})
names(cpgset) = paste0('chr',seq_along(cpgset))

show(cpgset);

noncpg = noncpgSitesFromCpGset(cpgset, 3);

show(noncpg);

isocpg = isocpgSitesFromCpGset(cpgset, 3);

show(isocpg);


}
\keyword{ bam }
\keyword{ Rbam }
