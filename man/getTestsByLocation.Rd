\name{getTestsByLocation}
\alias{getTestsByLocation}
\title{
Extract MWAS Test Results by Genomic Location
}
\description{
    Extract MWAS test results by genomic location.
}
\usage{
getTestsByLocation(x, chr, position)
}
\arguments{
    \item{x}{
        Parameter list or MWAS directory.
    }
    \item{chr}{
        Chromosome name or number.
    }
    \item{position}{
        Position(s) on the chromosome.
    }
}
\details{
    The function returns only records for locations that were tested.
    Other locations are ignored.
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
# Extract one location
getTestsByLocation(param, 1, 123321)

# Chromosome can be character
getTestsByLocation(param, "chr1", 123321)

# Extract multiple locations
# on the same chromosome
getTestsByLocation(param, 1, c(123123,123321))

# on different chromosomes
getTestsByLocation(param, c("chr1","chr2"), c(123123,123321))
}
}
