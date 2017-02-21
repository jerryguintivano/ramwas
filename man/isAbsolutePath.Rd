\name{isAbsolutePath}
\alias{isAbsolutePath}
\title{
    Check if Path is Absolute.
}
\description{
    Check whether a path is relative or absolute.
}
\usage{
isAbsolutePath(path)
}
\arguments{
    \item{path}{
        Path to be tested.
    }
}
\details{
    The function is designed to word with both Windows and Unix paths.
}
\value{
    \code{TRUE} if the path is absolute, \code{FALSE} otherwise.
}
\author{
    Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\note{
    This function improves upon the analog function
    in \code{R.utils} package.
    For instance, "~hi" is not an absolute path.
}
\seealso{
    See also \code{\link{makefullpath}}.
}
\examples{
isAbsolutePath( "C:/123" )   # TRUE
isAbsolutePath( "~123" )     # FALSE
isAbsolutePath( "~/123" )    # TRUE
isAbsolutePath( "/123" )     # TRUE
isAbsolutePath( "\\\\123" )    # TRUE
isAbsolutePath( "asd\\\\123" ) # FALSE
isAbsolutePath( "a\\\\123" )   # FALSE
}
