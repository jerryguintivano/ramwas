\name{findBestNpvs}
\alias{findBestNpvs}
\title{
    Quickly Find N Smallest P-values in a Long Vector
}
\description{
    Finding top, say, 100 p-values out of millions can be slow.
    This function does it much faster than the usual application of
    \code{order(pv)[1:N]}.
}
\usage{
findBestNpvs(pv, n)
}
\arguments{
    \item{pv}{
        Vector of p-values.
    }
    \item{n}{
        Number of best p-values to select.
    }
}
\details{
    The function is a faster analog of \code{sort(order(pv)[1:N])}
}
\value{
    Return a vector of positions of the smallest \code{N} p-values in \code{pv}.
}
\author{
    Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\seealso{
    See \code{\link[base]{order}}.
}
\examples{
pv = runif(1000)^10
n = 100

# Faster version
topSites1 = findBestNpvs(pv, n)

# Slow alternative
topSites2 = sort(order(pv)[1:n])

# The results must match
stopifnot(all( topSites1 == topSites2 ))
}
