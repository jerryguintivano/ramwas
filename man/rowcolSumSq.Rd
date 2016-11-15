\name{rowcolSumSq}
\alias{rowSumsSq}
\alias{colSumsSq}
\title{
    Form Row and Column Sums of Squares
}
\description{
    Form row and column sums of squares for numeric matrices.
    The functions are introduced as faster analogs of
    \code{rowSums(x^2)} and \code{colSums(x^2)} calls.
}
\usage{
rowSumsSq(x)
colSumsSq(x)
}
\arguments{
    \item{x}{
        Numeric matrix.
    }
}
\details{
    The function is implemented in C for better performance.
}
\value{
    Return a vector of sums of values in each row/column for matrix
    \code{x} (\code{rowSumsSq}/\code{colSumsSq}).
}
\author{
    Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\seealso{
    See \code{\link[base]{rowSums}} and \code{\link[base]{colSums}}
    for simple (not squared) row/column sums.
}
\examples{
x = matrix( 1:99, 9, 11)


# Calculate sums of squared elements in each row
rsum2 = rowSumsSq(x)

# Compare with alternative calculation
stopifnot( all.equal( rsum2, rowSums(x^2) ))


# Calculate sums of squared elements in each column
csum2 = colSumsSq(x)

# Compare with alternative calculation
stopifnot( all.equal( csum2, colSums(x^2) ))
}
