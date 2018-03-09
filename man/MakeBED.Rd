\name{madeBED}
\alias{madeBED}
\alias{madeBEDrange}
\alias{madeBEDgraph}
\alias{madeBEDgraphRange}
\title{
    Export MWAS results in BED format.
}
\description{
    Functions for exporting MWAS results in BED format files.
    
    Function \code{madeBED} saves MWAS findings in BED format for 
    all variables (CpGs), while \code{madeBEDrange} selects only variables on 
    a given chromosome between given locations.
    
    Functions \code{madeBEDgraph} and \code{madeBEDgraphRange} do the same,
    but create a file in BedGraph format.
}
\usage{
madeBED(x, filename)
madeBEDrange(x, filename, chr, start, end)
madeBEDgraph(x, filename)
madeBEDgraphRange(x, filename, chr, start, end)
}
\arguments{
    \item{x}{
        Name of MWAS directory (parameter \code{dirmwas}) or   
        list of RaMWAS parameters as described in the "RW6_param.Rmd" vignette.
        \cr
        Try: \code{vignette("RW6_param","ramwas")}.\cr
    }
    \item{filename}{
        Name of the BED file to create. If file exists, it's overwritten.
    }
    \item{chr}{
        Chromosome name or number.
    }
    \item{start}{
        Start position of the genomic region of interest.
    }
    \item{end}{
        End position of the genomic region of interest.
    }
}
\details{
    The function returns the MWAS results with locations.
}
\value{
    Returns a data.frame with BED file content:
    \item{chrom}{Chromosome}
    \item{chromStart}{Start position}
    \item{chromEnd}{End position}
    \item{name}{Empty name column. BED format only}
    \item{score}{p-value}
    
}
\author{
    Andrey A Shabalin \email{andrey.shabalin@gmail.com}
}
\seealso{
    See vignettes: \code{browseVignettes("ramwas")}.
}
\examples{
\dontrun{
# Extract BED file using parameter vector
madeBED(param, "file.bed")
madeBEDrange(param, "file.bed", 1, 123321, 223321)

# Extract BED file using directory name
madeBED("/data/myMWAS", "file.bed")
madeBEDrange("/data/myMWAS", "file.bed", 1, 123321, 223321)
}
}
