\name{parameterPreprocess}
\alias{parameterPreprocess}
\title{
    Preprocess Pipeline Parameter List.
}
\description{
    Fill in missing parameters with default values,
    read supporting data files,
    make relative directory path parameters absolute.
}
\usage{
parameterPreprocess(param)
}
\arguments{
    \item{param}{
        List with RaMWAS parameters.\cr
        For detailed description of all available parameters run:\cr
        \code{browseVignettes("ramwas")}.
    }
}
\details{
    A number of common preprocessing steps necessary for parameters of
    multiple pipeline parts are combined in this function.
    The actions include
    \itemize{
        \item Fill in default values for all missing parameters.
        \item Set \var{bamnames} parameter to the content
        \var{filebamlist} file (if \var{bamnames} was not set).
        \item Set \var{bam2sample} parameter to processed content of
        \var{filebam2sample} file (if \var{bam2sample} was not set).
        \item Set \var{covariates} parameter to the data frame from
        \var{filecovariates} file (if \var{covariates} was not set).
        \item Check parameters for consistency, i.e. that
        \var{modelcovariates} include only names of columns in
        \var{covariates}.
        \item Check that files \var{filecpgset} and
        \var{filenoncpgset} exist if the parameters are set.
    }
}
\value{
    Returns preprocessed list of parameters.
}
\author{
    Andrey A Shabalin \email{andrey.shabalin@gmail.com}
}
\note{
    This function is not intended to be run by the user.
}
\seealso{
    See vignettes: \code{browseVignettes("ramwas")}.
}
\examples{
param = ramwasParameters(
    dirproject = "."
)

param2 = parameterPreprocess(param)

print(param2)
}
