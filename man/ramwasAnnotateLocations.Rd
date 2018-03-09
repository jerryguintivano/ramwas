\name{ramwasAnnotateLocations}
\alias{ramwasAnnotateLocations}
\title{
    Extract Biomart Annotation for a Vector of Locations.
}
\description{
    Calls biomart annotation database
    for a vector of locations and
    assignes the tracks to the locations.
}
\usage{
ramwasAnnotateLocations(param, chr, pos)
}
\arguments{
    \item{param}{
        List of parameters as described in the "RW6_param.Rmd" vignette. \cr
        Try: \code{vignette("RW6_param","ramwas")}.
    }
    \item{chr}{
        A vector of chromosome names or numbers.
    }
    \item{pos}{
        A vector of genomic locations on the chromosomes.
    }
}
\details{
    This function is used internally by RaMWAS annotation step.
}
\value{
    An annotation table, on line per supplied location.
}
\author{
    Andrey A Shabalin \email{andrey.shabalin@gmail.com}
}
\seealso{
    See vignettes: \code{browseVignettes("ramwas")}
}
\examples{
bihost = "grch37.ensembl.org"
bimart = "ENSEMBL_MART_ENSEMBL"
bidataset = "hsapiens_gene_ensembl"
biattributes = c("hgnc_symbol", "entrezgene", "strand")
bifilters = list(with_hgnc_trans_name = TRUE)
biflank = 0

param = ramwasParameters(
    bihost = bihost,
    bimart = bimart,
    bidataset = bidataset,
    biattributes = biattributes,
    bifilters = bifilters,
    biflank = biflank)

# Test a location
chr = "chr1"
pos =  15975530
\dontrun{
ramwasAnnotateLocations(param, chr, pos)
}
}
