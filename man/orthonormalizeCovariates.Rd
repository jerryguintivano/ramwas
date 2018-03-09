\name{orthonormalizeCovariates}
\alias{orthonormalizeCovariates}
\title{
    Orthonormalize Covariates
}
\description{
    Takes a matrix of data frame with covariates,
    adds a constant covariate (optional), and
    orthonormalizes the set.
}
\usage{
orthonormalizeCovariates(cvrt, modelhasconstant)
}
\arguments{
    \item{cvrt}{
        A matrix or data frame with covariates (one column per covariate).
    }
    \item{modelhasconstant}{
        Set to \code{TRUE} to add a constant covariate
        into the set before normalization.
    }
}
\details{
    Factor variables are split into dummy variables before orthonormalization.
    \cr
    The operation is performed via QR decomposition (\link[base]{qr}).
}
\value{
    Returns a matrix with orthogonal columns with unit length,
    whose columns spans the same space as the covariates plus a constant
    (if \code{modelhasconstant} is \code{TRUE}).
}
\author{
    Andrey A Shabalin \email{andrey.shabalin@gmail.com}
}
\note{
    This function is used in several parts of the pipeline.
}
\examples{
# Sample matrix of covariates
covariates = data.frame(a = 1:12, b = 12:1)

# Orthonormalizing Covariates
cvrtqr = orthonormalizeCovariates(covariates, modelhasconstant = TRUE)

# Checking the results (round to ignore rounding errors)
print( round(crossprod(cvrtqr),15) )

# Stop if not orthonormal
stopifnot(all.equal( crossprod(cvrtqr), diag(ncol(cvrtqr)) ))

# Example with a factor variable
groups = data.frame(gr = c("a","a","a","b","b","b","c","c","c"))
orthonormalizeCovariates(groups)
}
