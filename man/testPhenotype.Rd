\name{testPhenotype}
\alias{testPhenotype}
\title{
    Test the Phenotype of Interest for Association with Methylation Coverage.
}
\description{
    An internal, function for fast association testing.
    It tests the phenotype of interest for
    association with methylation coverage (columns of the data parameter).
}
\usage{
testPhenotype(phenotype, data1, cvrtqr)
}
\arguments{
    \item{phenotype}{
        Vector with phenotype. Can be numerical, character, or factor vector.
    }
    \item{data1}{
        Matrix with data (normalized coverage), one variable (CpG) per column.
    }
    \item{cvrtqr}{
        Orthonormalized covariates, one covariate per column.
        See \code{\link{orthonormalizeCovariates}}.
    }
}
\details{
    The testing is performed using matrix operations and C/C++ code,
    emplying an approach similar to that in \pkg{MatrixEQTL}.
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
    \item{nVarTested}{First number of degrees of freedom for the F-test.
    Equal to the number of factor levels reduced by 1}
    \item{dfFull}{Second number of degrees of freedom for the F-test.}
}
\author{
    Andrey A Shabalin \email{andrey.shabalin@gmail.com}
}
\note{
    This function is used in several parts of the pipeline.
}
\seealso{
    See vignettes: \code{browseVignettes("ramwas")}.

    Also check \code{\link{orthonormalizeCovariates}}.
}
\examples{
### Generate data inputs
# Random data matrix with signal in the first column
data = matrix(runif(30*5), nrow = 30, ncol = 5)
data[,1] = data[,1] + rep(0:2, each = 10)

# Two random covariates
cvrt = matrix(runif(2*30), nrow = 30, ncol = 2)
cvrtqr = orthonormalizeCovariates(cvrt)



### First, illustrate with numerical phenotype
# Numerical, 3 value phenotype
phenotype = rep(1:3, each = 10)

# Test for association
output = testPhenotype(phenotype, data, cvrtqr)

# Show the results
print(output)

# Comparing with standard R code for the first variable
summary(lm( data[,1] ~ phenotype + cvrt ))



### First, illustrate with numerical phenotype
# Categorical, 3 group phenotype
phenotype = rep(c("Normal", "Sick", "Dead"), each = 10)

# Test for association
output = testPhenotype(phenotype, data, cvrtqr)

# Show the results
print(output)

# Comparing with standard R code for the first variable
anova(lm( data[,1] ~ cvrt + phenotype ))
}
