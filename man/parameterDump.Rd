\name{parameterDump}
\alias{parameterDump}
\title{
    Save Parameters in a Text File
}
\description{
    Saves parameters in a text file,
    prioritizing those listed in \code{toplines}.
}
\usage{
parameterDump(dir, param, toplines = NULL)
}
\arguments{
    \item{dir}{
        Directory to save the parameters to.
        The file is named \code{"UsedSettings.txt"}.
    }
    \item{param}{
        A list with RaMWAS parameters. Or any list in general.
    }
    \item{toplines}{
        Names of the elements in \code{param} to save first (top of the file).
    }
}
\details{
    This function is used internally by multiple RaMWAS
    functions to record parameters used to run the analysis.
}
\value{
    The function creates a file and returns nothing.
}
\author{
    Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\note{
    This function is not intended to be run by the user.
}
\seealso{
    See vignettes: \code{browseVignettes("ramwas")}
}
\examples{
param = ramwasParameters(
    number = 123123,
    integer = 312L,
    textline = "Hi there",
    characterVector = c("Hi","Hi again","Bye"),
    dataframe = data.frame(a = 1:12, b = 12:1)
    )

thedir = tempdir()

parameterDump(thedir, param, c("integer","characterVector"))

cat( readLines( paste0(thedir,"/UsedSettings.txt") ), sep = "\n")

file.remove( paste0(thedir,"/UsedSettings.txt") )
}
