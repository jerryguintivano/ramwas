\name{plotFragmentSizeDistributionEstimate}
\alias{plotFragmentSizeDistributionEstimate}
\alias{estimateFragmentSizeDistribution}
\title{
    Estimate and plot Fragment Size Distribution.
}
\description{
    RaMWAS functions for estimation and plotting of
    the fragment size distribution.
}
\usage{
estimateFragmentSizeDistribution(frdata, seqLength)
plotFragmentSizeDistributionEstimate(      
        frdata,
        estimate,
        col1 = "blue",
        col2 = "red")
}
\arguments{
    \item{frdata}{
        Distribution of distances from the starts of isolated reads
        to the respective CpGs.
    }
    \item{seqLength}{
        The length of sequenced part of the fragments.\cr
        The fragments are assument to not be smaller than \code{seqLength}.
    }
    \item{estimate}{
        Fragment size distribution estimate.
    }
    \item{col1}{
        Color of \code{frdata} points.
    }
    \item{col2}{
        Color of \code{estimate} curve.
    }
}
\value{
    The function \code{estimateFragmentSizeDistribution}
    returns the estimate of the fragment size distribution.
}
\author{
    Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\note{
    If the length of \code{frdata} is equal to \code{seqLength},
    the fragments are assumed to all be of length \code{seqLength}.
}
\seealso{
    See vignettes: \code{browseVignettes("ramwas")}.
}
\examples{
# Simulate data
x = 0:250
truemean = 1 - pnorm(x, mean = 150, sd = 50)
frdata = rpois(n = length(x), lambda = truemean*300)

# Estimate fragment size distribution
estimate = estimateFragmentSizeDistribution(frdata, seqLength = 50)

# Plot fragment size distribution estimate
plotFragmentSizeDistributionEstimate(frdata, estimate)
}
