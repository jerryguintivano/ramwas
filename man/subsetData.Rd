\name{subsetData}
\alias{subsetCoverageDirByLocation}
\title{
    Subset a data matrix and locations
}
\description{
    Subset a data (coverage) matrix and corresponding matrix of locations
    to a specified set of locations.
}
\usage{
subsetCoverageDirByLocation(x, chr, start, targetdir)
}
\arguments{
    \item{x}{
        Name of data (coverage) directory or   
        list of RaMWAS parameters as described in the "RW6_param.Rmd" vignette.
        \cr
        Try: \code{vignette("RW6_param","ramwas")}.
    }
    \item{chr}{
        Vector of chromosome names or numbers.
    }
    \item{start}{
        Start positions of the CpGs of interest.
    }
    \item{targetdir}{
        Directory name for the new (subset) data matrix and locations.
    }
}
\value{
    The function returns nothing.
}
\author{
    Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\seealso{
    See vignettes: \code{browseVignettes("ramwas")}.
}
\examples{

x = "/data/myCoverageMatrix"
chr = c("chr1", "chr2", "chr3")
start = c(12345, 123, 12)
targetdir = "/data/subsetCoverageMatrix"

\dontrun{
subsetCoverageDirByLocation(x, chr, start, targetdir)
}
}
