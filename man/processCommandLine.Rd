\name{processCommandLine}
\alias{processCommandLine}
\title{
    Scan Parameters From Command Line
}
\description{
    The pipeline parameters can be provided via command line.

    For example: \cr
    \code{R pipeline.r dirproject="/project" maxrepeats=0 modeloutcome="Age"}

    Each command line argument is treated as an R statement.

    All variables defined this way are collected in a list which is returned.
}
\usage{
processCommandLine(.arg = NULL)
}
\arguments{
    \item{.arg}{
        Vector of command line parameters. Obtained from
        \code{\link[base]{commandArgs}} if omitted.
    }
}
\details{
    If a command line argument defines variable \code{"fileparam"},
    it is assumed to be a filename, and the file with this name
    is scanned for extra pipeline parameters,
    as by \code{\link{parametersFromFile}}.
}
\note{
    Variables with names starting with period (.) are ignored.
}
\value{
    Returns the list with all the variables set by
    the statement in the command line.
}
\author{
    Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\seealso{
    See vignettes: \code{browseVignettes("ramwas")}.
}
\examples{

filename = tempfile()

# Assume command line with two components:
# dirproject="."
# modelcovariates=c("Age","Sex")

arg = c(
    "dirproject = \".\"",
    "modelcovariates = c(\"Age\",\"Sex\")")

# Process the command line
param = processCommandLine(arg)

# Show the list
print(param)
}
