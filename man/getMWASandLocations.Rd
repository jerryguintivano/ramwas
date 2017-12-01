\name{getMWASandLocations}
\alias{getMWASandLocations}
\title{
Extract MWAS Test Results for All Genomic Location
}
\description{
    Extract MWAS test results for all genomic location.
}
\usage{
getMWASandLocations(x)
}
\arguments{
    \item{x}{
        Parameter list or the MWAS directory.
    }
}
\details{
    The function returns the MWAS results with locations.
}
\value{
    Returns a data frame with
    \item{chr}{Chromosome}
    \item{position}{position}
    \item{cor}{coverage - phenotype correlation}
    \item{t.test}{t-statistic}
    \item{p.value}{p-value}
    \item{q.value}{q-value (FDR)}

    If the outcome variable was categorical,
    columns \code{R.squared} and \code{F-test} are present
    in place of \code{cor} and \code{t.test}.
}
\author{
    Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\seealso{
    See vignettes: \code{browseVignettes("ramwas")}.
}
\examples{
\dontrun{
# Extract MWAS using parameter vector
getMWASandLocations(param)

# Extract MWAS using directory name
getMWASandLocations("/data/myMWAS")
}
}
