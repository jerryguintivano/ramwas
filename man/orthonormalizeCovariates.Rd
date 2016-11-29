\name{orthonormalizeCovariates}
\alias{orthonormalizeCovariates}
\title{
    Orthonormalize Covariates
}
\description{
    Takes a set of covariates, add a constant covariate, and
    orthonormalizes the set.
}
\usage{
orthonormalizeCovariates(cvrt)
}
\arguments{
    \item{cvrt}{
        A matrix or data frame with covariates.
    }
}
\details{
    Factor variables are split into dummy variables before orthonormalization.

    The operation is performed via QR decomposition (\link[base]{qr}).
}
\value{
    A matrix with orthogonal columns with unit length,
    whose columns spans the same space as the covariates plus a constant.
}
\author{
    Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\note{
    This function is used in several parts of the pipeline.
}
\examples{
# Sample matrix of covariates
covariates = data.frame(a = 1:12, b = 12:1)

# Orthonormalizing Covariates
cvrtqr = orthonormalizeCovariates(covariates)

# Checking the results (round to ignore rounding errors)
print( round(crossprod(cvrtqr),15) )

# Stop if not orthonormal
stopifnot(all.equal( crossprod(cvrtqr), diag(ncol(cvrtqr)) ))
}
