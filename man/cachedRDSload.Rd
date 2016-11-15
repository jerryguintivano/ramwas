\name{cachedRDSload}
\alias{cachedRDSload}
\title{Cached Loading of RDS Files}
\description{
Loads an .rds file \code{rdsfilename} using \code{\link[base]{readRDS}}
and returns the loaded object.
The object is also saved in a cache so that repeated calls of the
function with the same filename return the same object instanteneously.
}
\usage{
cachedRDSload(rdsfilename)
}
\arguments{
    \item{rdsfilename}{
        Name of the RDS file.
    }
}
\details{
    The cached object is stored in a private package environment.
}
\value{
    Returns the object loaded with \code{\link[base]{readRDS}} from
    \code{rdsfilename} at this or a previous call of the function.
}
\author{
    Andrey Shabalin \email{ashabalin@vcu.edu}
}
\examples{
### Change filename to hg19 CpGset

filename = system.file("extdata", "qc_sample.rds", package = "ramwas")

time1 = system.time( {obj1 = cachedRDSload(filename)} )
time2 = system.time( {obj1 = cachedRDSload(filename)} )

cat("First  loading time:",time1[3],"seconds","\n")
cat("Second loading time:",time2[3],"seconds","\n")
}
\keyword{ bam }
\keyword{ Rbam }
