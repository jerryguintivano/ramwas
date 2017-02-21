\name{pipeline}
\alias{pipeline}
\alias{ramwas1scanBams}
\alias{pipelineProcessBam}
\alias{ramwas2collectqc}
\alias{ramwas3normalizedCoverage}
\alias{ramwas4PCA}
\alias{ramwas5MWAS}
\alias{ramwas6annotateTopFindings}
\alias{ramwas7ArunMWASes}
\alias{ramwas7BrunElasticNet}
\alias{ramwas7CplotByNCpGs}
\alias{ramwas7riskScoreCV}
\alias{ramwasSNPs}
\title{
    RaMWAS: High Level Pipeline Functions
}
\description{
    These functions provide a simple way to run all steps of RaMWAS pipeline.
}
\usage{
ramwas1scanBams(param)
pipelineProcessBam(bamname, param)
ramwas2collectqc(param)
ramwas3normalizedCoverage(param)
ramwas4PCA(param)
ramwas5MWAS(param)
ramwas6annotateTopFindings(param)
ramwas7ArunMWASes(param)
ramwas7BrunElasticNet(param)
ramwas7CplotByNCpGs(param)
ramwas7riskScoreCV(param)
ramwasSNPs(param)
}
\arguments{
    \item{param}{
        List with RaMWAS parameters. \cr
        For detailed description of all available parameters run:\cr
        \code{browseVignettes("ramwas")}.
        }
    \item{bamname}{
        Name of the BAM file to process.
        Can be absolute or relative to \code{dirbam} parameter
        (in \code{param} list).
    }
}
\details{
    See vignettes for details: \code{browseVignettes("ramwas")}.
}
\value{
    Function \code{pipelineProcessBam} returns "OK. <bamname>"
    if no error occurred. Otherwise, returns text with error.
    Other functions return nothing.
}
\author{
    Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\seealso{
    See vignettes: \code{browseVignettes("ramwas")}.
}
\examples{
param = ramwasParameters(
    dirbam = "/project/bams",
    dirproject = "/project",
    filebamlist = "000_list_of_files.txt",
    scoretag = "AS",
    minscore = 100,
    cputhreads = 4,
    filecpgset    = "/RaMWAS/hg19_1kG_MAF_0.01_chr1-22_bowtie2_75bp.rds",
    filenoncpgset = "/RaMWAS/hg19_1kG_MAF_0.01_chr1-22_bowtie2_75bp_nonCpG.rds",
    maxrepeats = 3,
    maxfragmentsize = 250,
    minfragmentsize = 75,
    filebam2sample = "000_list_of_files.txt",
    filecovariates = "Covariates.txt",
    modelcovariates = c("Age","Sex"),
    modeloutcome = "CellType",
    modelPCs = 1,
    cvnfolds = 10,
    mmncpgs = 1000,
    mmalpha = 0
)

\dontrun{
ramwas1scanBams(param)
ramwas2collectqc(param)
ramwas3normalizedCoverage(param)
ramwas4PCA(param)
ramwas5MWAS(param)
ramwas6annotateTopFindings(param)
ramwas7riskScoreCV(param)}
}
