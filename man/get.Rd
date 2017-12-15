\name{get}
\alias{getLocations}
\alias{getMWAS}
\alias{getMWASandLocations}
\alias{getTestsByLocation}
\alias{getDataByLocation}
\title{
    Functions for Access to Data, MWAS Results, and Location Information
}
\description{
    Functions for access to data, MWAS results, and location information.
    
    Function \code{getLocations} obtains the location information for
    all variables (CpGs).
    
    Function \code{getMWASandLocations} obtains both MWAS results and
    location information in a single data frame.
    
    Functions \code{getDataByLocation} and \code{getTestsByLocation} 
    return the data (coverage) and MWAS results for the selected set
    of variables (CpGs).
}
\usage{
getLocations(x)
getMWAS(x)
getMWASandLocations(x)
getTestsByLocation(x, chr, position)
getDataByLocation(x, chr, position)
}
\arguments{
    \item{x}{
        Name of directory or   
        list of RaMWAS parameters as described in the "RW6_param.Rmd" vignette. \cr
        Try: \code{vignette("RW6_param","ramwas")}.\cr

        If a directory name is provided, it must point to
        \itemize{
            \item Data (coverage) directory (parameter \code{dircoveragenorm})
                for \code{getDataByLocation} and \code{getLocations}
            \item MWAS directory (parameter \code{dirmwas}) 
                for \code{getMWAS}, \code{getMWASandLocations},
                and \code{getTestsByLocation}
        }
    }
    \item{chr}{
        vector of chromosome names or numbers.
    }
    \item{position}{
        Start positions on the respective chromosome(s).
    }
}
\details{
    The function returns the MWAS results with locations.
}
\value{
    Function \code{getLocations} returns a data frame with
    \item{chr}{Chromosome}
    \item{start}{Start position}
    \item{end}{End position}
    
    Function \code{getMWAS} returns a data frame with
    \item{cor}{coverage - phenotype correlation}
    \item{t.test}{t-statistic}
    \item{p.value}{p-value}
    \item{q.value}{q-value (FDR)}

    If the outcome variable was categorical,
    columns \code{cor} and \code{t.test} are replaced with
    \code{R.squared} and \code{F-test}

    Functions \code{getMWASandLocations} and \code{getTestsByLocation}
    return a data frame with elements of output
    of both \code{getLocations} and \code{getMWAS}

    
    Function \code{getDataByLocation} returns a list with
    \item{locations}{Chromosomal location information for located variables}
    \item{matrix}{Data (coverage) matrix for the selected locations}
}
\author{
    Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\seealso{
    See vignettes: \code{browseVignettes("ramwas")}.
}
\examples{
\dontrun{
# Extract locations using parameter vector
getLocations(param)

# Extract locations using directory name
getLocations("/data/myMWAS")


# Extract MWAS using parameter vector
getMWAS(param)

# Extract MWAS using directory name
getMWAS("/data/myMWAS")


# Extract MWAS using parameter vector
getMWASandLocations(param)

# Extract MWAS using directory name
getMWASandLocations("/data/myMWAS")


# Extract one location
getTestsByLocation(param, 1, 123321)

# Chromosome can be character
getTestsByLocation(param, "chr1", 123321)

# Extract multiple locations
# on the same chromosome
getTestsByLocation(param, 1, c(123123,123321))

# on different chromosomes
getTestsByLocation(param, c("chr1","chr2"), c(123123,123321))


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
