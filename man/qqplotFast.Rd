\name{qqPlotFast}
\alias{qqPlotFast}
\title{
    Fast QQ-plot for Large Number of P-values
}
\description{
    Plots a QQ-plot with a confidence band and
    an estimate of inflation factor lambda.
    Works quickly even for tens of millions of p-values.
}
\usage{
qqPlotFast(pvalues, ntests = NULL, ci.level = 0.05)
}
\arguments{
    \item{pvalues}{
        Vector of p-values, preferably sorted.
    }
    \item{ntests}{
        If only significant p-values are provided,
        the total number of tests performed.
        By default equal to the length of \code{pvalues}.
    }
    \item{ci.level}{
        Significance level of the confidence band.
    }
}
\details{
    The plot has no title. To add a title use \code{\link[graphics]{title}}.
}
\value{
    This function plots a QQ-plot and returns nothing.
}
\author{
    Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\note{
    The function works faster if the p-values are sorted.
}
\seealso{
    See vignettes: \code{browseVignettes("ramwas")}.
}
\examples{
# Million p-values
n = 1e6

# Null p-values
pv = runif(n)

# QQ-plot should be nearly diagonal
qqPlotFast(pv)
}
