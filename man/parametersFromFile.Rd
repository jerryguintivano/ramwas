\name{parametersFromFile}
\alias{parametersFromFile}
\title{
    Scan Parameters From a R Code File
}
\description{
    The pipeline parameters can be stored in a simple file,
    formatted as R code.
    The \code{parametersFromFile} function transforms them into
    a parameter list used by RaMWAS steps.
}
\usage{
parametersFromFile(.parameterfile)
}
\arguments{
    \item{.parameterfile}{
        Name of the file with the parameters set as R variables.
        See the example below.
    }
}
\details{
    Variables with names starting with period (.) are ignored.
}
\note{
    The file \code{.parameterfile} is executed as R code,
    so use only trusted parameter files.
}
\value{
    Returns the list with all the variables set in the file.
}
\author{
    Andrey A Shabalin \email{andrey.shabalin@gmail.com}
}
\seealso{
    See vignettes: \code{browseVignettes("ramwas")}.
}
\examples{

filename = tempfile()

# Create a file with lines
# dirproject = "."
# modelcovariates = c("Age","Sex")

writeLines(
    con = filename,
    text = c(
        "dirproject = \".\"",
        "modelcovariates = c(\"Age\",\"Sex\")")
    )

# Scan the file into a list
param = parametersFromFile(filename)

# Show the list
print(param)

file.remove(filename)
}
