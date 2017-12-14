\name{plotPC}
\alias{plotPCvalues}
\alias{plotPCvectors}
\title{
    Plot Principal component (PC) Values (variation explained) and 
    PC vectors (loadings)
}
\description{
    The function \code{plotPCvalues} plots PC values (variation explained).\cr
    The function \code{plotPCvectors} plots PC vectors (loadings).
}
\usage{
plotPCvalues(values, n = 40)
plotPCvectors(e, i)
}
\arguments{
    \item{values}{
        Vector of PC values.
    }
    \item{n}{
        Number of top PCs to plot.
    }
    \item{e}{
        Eigenvalue decomposition (spectral decomposition) as returned by 
        \code{\link[base]{eigen}}.
    }
    \item{i}{
        Indicates loadings of which PC to plot.
    }
}
\value{
    This function creates a PC plot and returns nothing (\code{NULL}).
}
\author{
    Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\seealso{
    See vignettes: \code{browseVignettes("ramwas")}.
}
\examples{
# Sample data
# for 1000 observations and 10 samples
m = 1000
n = 10
data = matrix(rnorm(n*m), nrow = m)

# Covariance and eigenvalue decomposition
covmat = crossprod(data)
e = eigen(covmat)

# Plot PC values
plotPCvalues(e$values)

# Plot PC vectors
plotPCvectors(e, 1)
plotPCvectors(e, 2)
}
