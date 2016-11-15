\name{testPhenotype}
\alias{testPhenotype}
\title{
    Test the Phenotype of Interest for Association with Methylation Coverage.
}
\description{
    An internal, low-level function for fast association testing.
    It tests the phenotype of interest for
    association with methylation coverage.
}
\usage{
testPhenotype(phenotype, data, cvrtqr)
}
\arguments{
    \item{phenotype}{
        Vector with phenotype. Can be numerical or character/factor variable.
    }
    \item{data}{
        Matrix with normalized coverage, one CpG per column.
    }
    \item{cvrtqr}{
        Orthonormalized covariates (transposed).
        See \code{\link{orthonormalizeCovariates}}.
    }
}
\details{
    The testing is performed in a fast way,
    using an approach similar to that in \pkg{MatrixEQTL}.
}
\value{
    If the phenotype is numerical, the output is a list with
    \item{correlation}{Correlations between residualized phenotype and
    data columns.}
    \item{tstat}{Corresponding T-statistics}
    \item{pvalue}{Corresponding P-values}
    \item{nVarTested}{Always 1}
    \item{dfFull}{Number of degrees of freedom of the T-test}

    If the phenotype is a factor (or character)
    \item{Rsquared}{R-squared for the residualized ANOVA F-test.}
    \item{Fstat}{Corresponding F-test}
    \item{pvalue}{Corresponding P-values}
    \item{nVarTested}{First number of degrees of freedom of the F-test.
    Equal to the number of factor levels reduced by 1}
    \item{dfFull}{Second number of degrees of freedom of the F-test.}
}
\author{
    Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\note{
    This function is used in several parts of the pipeline.
}
\seealso{
    See vignettes: \code{browseVignettes("ramwas")}.

    Also check \code{\link{orthonormalizeCovariates}}.
}
\examples{
# Random data with signal in the first column
data = matrix( runif(30*3), 30, 3)
data[,1] = data[,1] + rep(0:2, each = 10)

# Random covariate
cvrtqr = orthonormalizeCovariates(matrix(runif(30),ncol=1))


# Categorical, 3 group phenotype
phenotype = rep(c("Normal","Sick","Dead"), each = 10)

# Test for association
output = testPhenotype(phenotype, data, t(cvrtqr))

print(output)

# Numerical, 3 value phenotype
phenotype = rep(1:3, each = 10)

# Test for association
output = testPhenotype(phenotype, data, t(cvrtqr))

print(output)
}
