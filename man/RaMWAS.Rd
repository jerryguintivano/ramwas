\name{ramwas-package}
\alias{ramwas-package}
\alias{ramwas}
\docType{package}
\title{
	Fast methylome-wide association study pipeline for enrichment platforms
}
\description{
Ramwas provides a full set of tools for 
 methylome-wide association studies for enrichment technologies.
 The steps of the pipeline include scanning aligned reads from BAM files,
 calculation of multiple quality control measures,
 creation of coverage for each sample,
 principal component analysis for capturing batch effects and accounting for 
 cell composition differences, and, finally, association analysis
 with respect to phenotypes of interest correcting for PCs and known covariates.
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
\references{
	.%The package website: \url{http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/}
}
\keyword{ package }
\keyword{ ramwas }
\keyword{ RaMWAS }
\seealso{
	Run \code{browseVignettes('ramwas')} for the vignettes.
}
