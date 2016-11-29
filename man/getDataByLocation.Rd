\name{getDataByLocation}
\alias{getDataByLocation}
\title{
Extract Data (Coverage) Matrix Columns by Genomic Location
}
\description{
    Extract data (coverage) matrix columns by genomic location.
}
\usage{
getDataByLocation(x, chr, position)
}
\arguments{
    \item{x}{
        Parameter list or directory with coverage matrix.
    }
    \item{chr}{
        Chromosome name or number.
    }
    \item{position}{
        Position(s) on the chromosome.
    }
}
\details{
    The function returns only data for locations that were tested.
    Other locations are ignored.
}
\value{
    Returns a list with 3 entries
    \item{chr}{Chromosome}
    \item{position}{position}
    \item{mat}{Data matrix for selected genomic locations}
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
getDataByLocation(param, 1, 123321)

# Chromosome can be character
getDataByLocation(param, "chr1", 123321)

# Extract multiple locations
# on the same chromosome
getDataByLocation(param, 1, c(123123,123321))

# on different chromosomes
getDataByLocation(param, c("chr1","chr2"), c(123123,123321))
}
}
