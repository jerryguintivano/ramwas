\name{toycpgset}
\alias{toycpgset}
\docType{data}
\title{
	Sample set of CpGs in RAMWAS format.
}
\description{
	Ramwas works with coverage calculated over a set of CpGs.
	The set of CpGs is stored as a \code{\link[base]{list}} with one sorted vector of CpG locations per chromosome.
	
	The sample CpG set is defined to have 3-4 CpGs on chromosomes chr1, chr2, and chr3:
\preformatted{
toycpgset = list( chr1 = c(12L, 57L, 123L),
                  chr2 = c(45L, 95L, 99L, 111L),
                  chr3 = c(22L, 40L, 199L, 211L) );
}
}

\keyword{datasets}
