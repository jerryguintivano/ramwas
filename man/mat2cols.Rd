\name{mat2cols}
\alias{mat2cols}
\title{
    Split a Matrix into Column Vectors
}
\description{
    Internal function for splitting a matrix into column vectors.
}
\usage{
mat2cols(x)
}
\arguments{
    \item{x}{
        A matrix.
    }
}
\value{
    List of matrix columns.
}
\author{
    Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\seealso{
    See vignettes: \code{browseVignettes("ramwas")}.
}
\examples{
# Sample data
data = matrix(1:12, nrow = 3)

# Split it
mat2cols(data)
}
