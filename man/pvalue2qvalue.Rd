\name{pvalue2qvalue}
\alias{pvalue2qvalue}
\title{
    Calculate Benjamini-Hochberg q-values
}
\description{
    Calculate Benjamini-Hochberg q-values for a vector of p-values.
}
\usage{
pvalue2qvalue(pv, n = length(pv))
}
\arguments{
    \item{pv}{
        Vector of p-values.
    }
    \item{n}{
        If \code{pv} has only top p-values from a bigger set,
        \code{n} should indicate the number of tests performed.
    }
}
\details{
    The q-values can be slightly conservative compared to
    other popular q-value calculation methods.
}
\note{
    The function runs faster if the vector \code{pv} is sorted.
}
\value{
    Return a vector of q-values matching p-values in \code{pv}.
}
\author{
    Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\examples{
pv = runif(20)^2

qv = pvalue2qvalue(pv)
}
