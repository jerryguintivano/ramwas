\name{plotROC}
\alias{plotROC}
\title{
    Plot ROC Curve.
}
\description{
    Plot ROC Curve for a binary outcome and its predictor.
}
\usage{
plotROC(outcome, forecast)
}
\arguments{
    \item{outcome}{
        Outcome vector taking two distinct values
    }
    \item{forecast}{
        Numeric predictor vector.
        Larger values of \code{forecast} should correspond to
        more likely chance of the larger outcome.
    }
}
\details{
    The plot has no title. To add a title use \code{\link[graphics]{title}}.
}
\value{
    This function plots a ROC curve and returns nothing.
}
\author{
    Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\note{
    The input vectors must be of equal length and have matching values.
}
\seealso{
    See vignettes: \code{browseVignettes("ramwas")}.
}
\examples{
# Generate binary outcome and continuous predictor
# for 1000 observations
n = 1000;
outcome = sample(0:1, size = n, replace = TRUE)
forecast = outcome + rnorm(n);

# Make ROC plot
plotROC(outcome, forecast)
}
