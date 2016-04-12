\name{parameters}
\alias{parametersFromFile}
\title{Scan parameter file}
\description{
Scan parameter file as an R code (using \code{\link{source}}) and record all variables in a named list.
}
\usage{
parametersFromFile(.parameterfile)
}
\arguments{
	\item{.parameterfile}{
		Name of the R-formatted file with all the parameters.
	}
}
\details{
	.
}
\value{
	Returns a list with all the variables set in the \code{.parameterfile} file.
	Variables with names starting with period (.) are ignored.
}
\references{
	.
}
\seealso{
	.
}
\note{
	.
}
\author{
	Andrey Shabalin \email{ashabalin@vcu.edu}
}
\examples{
parameterfile = system.file("extdata", "param_file.txt", package = "ramwas")
param = parametersFromFile(parameterfile)

show(param)
}
\keyword{ bam }
\keyword{ Rbam }
