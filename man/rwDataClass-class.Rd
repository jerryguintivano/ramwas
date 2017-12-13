\name{rwDataClass-class}
\Rdversion{1.1}
\docType{class}
\alias{rwDataClass-class}
\alias{rwDataClass}

\title{Class for Accessing Data (Coverage) Matrix}
\description{
    This class is a wrapper for accessing the data (coverage) matrix.
    It automatically subsets the samples to those listed in the covariates.
    Data access function imputes missing values and 
    can residualize the variables.
}
\section{Extends}{
    \code{rwDataClass} is a reference classes 
    (see \code{\linkS4class{envRefClass}}).
}
\author{
    Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\seealso{
    See vignettes: \code{browseVignettes("ramwas")}
}
\keyword{classes}
\section{Fields}{
    \describe{
        \item{\code{fmdata}:}{
            Filematrix object for the data matrix. \cr
            Not intended to be accessed directly.
        }
        \item{\code{samplenames}:}{
            Vector of sample names.
        }
        \item{\code{nsamples}:}{
            Number of samples
        }
        \item{\code{ncpgs}:}{
            Number of variables (CpG sites) in the data matrix.
        }
        \item{\code{ndatarows}:}{
            Number of variables in the data matrix 
            (may be bigger than the number of samples).
        }
        \item{\code{rowsubset}:}{
            Indices of samples in the data matrix.
        }
        \item{\code{cvrtqr}:}{
            Matrix of orthonormalized covariates.
        }
    }
}
\section{Methods}{
    \describe{
        \item{\code{initialize(param = NULL, getPCs = TRUE, lockfile = NULL)}:}{
            Create the data access class.\cr
            'param' should contain the RaMWAS parameter vector.\cr
            'getPCs' indicates if the covariate set should include
            Principal components.\cr
            'lockfile' is the 'lockfile' parameter used in accessing the 
            data filematrix.
        }
        \item{\code{open(param = NULL, getPCs = TRUE, lockfile = NULL)}:}{
            The same as 'initialize' method, but for already created object.
        }
        \item{\code{close()}:}{
            Clears the object. Closes the filematrix.
        }
        \item{\code{getDataRez(colset, resid = TRUE)}:}{
            Extracts data for variables indexed by 'colset'.\cr
            The data is residualized unless \code{resid = FALSE}.
        }
    }
}
\examples{
# Create an empty rwDataClass
data = new("rwDataClass")

\dontrun{
# Connect to the data
data$open(param)

# Create a rwDataClass and connect to the data

data = new("rwDataClass", param = param)

# close the object
data$close()
}
}
