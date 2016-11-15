\name{ramwas-package}
\alias{ramwas-package}
\alias{ramwas}
\docType{package}
\title{
    Fast Methylome-wide Association Study Pipeline for Enrichment Platforms
}
\description{
RaMWAS provides a complete toolset for
    methylome-wide association studies (MWAS).
    It is specifically designed for data from
    enrichment based methylation assays,
    but can be applied to other methylomic data as well.
    The analysis pipeline includes seven steps:
    (1) scanning aligned reads from BAM files,
    (2) calculation of quality control measures,
    (3) creation of methylation score (coverage) matrix,
    (4) principal component analysis for capturing batch effects and
    detection of outliers,
    (5) association analysis with respect to phenotypes of interest
    while correctingfor top PCs and known covariates,
    (6) annotation of significant findings, and
    (7) multi-marker analysis (methylation risk score) using elastic net.
}
\details{
    \tabular{ll}{
        Package: \tab ramwas\cr
        Type: \tab Package\cr
        License: \tab LGPL-3 \cr
        LazyLoad: \tab yes\cr
        Depends: \tab methods\cr
    }
}
\author{
    Andrey Shabalin \email{ashabalin@vcu.edu}

    Maintainer: Andrey Shabalin <ashabalin@vcu.edu>
}
\keyword{ package }
\keyword{ ramwas }
\keyword{ RaMWAS }
\seealso{
    See vignettes: \code{browseVignettes("ramwas")}.
}
