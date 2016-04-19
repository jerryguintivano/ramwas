\name{pipeline}
\alias{pipelineProcessBam}
\alias{ramwas1scanBams}
\title{Cached loading of an .rds file.}
\description{
These functions provide a convenient, high level, way to run all steps of RaMWAS pipeline.
}
\usage{
pipelineProcessBam(bamname, param)

ramwas1scanBams(param)
}
\arguments{
	\item{param}{
		List of parameters for the RaMWAS pipeline.
	}
	\item{bamname}{
		File name for the BAM file. Can include relative or absolute path to the file.
	}
}
\details{
	\code{pipelineProcessBam} scans a single BAM file and records the read start locations and QC information.
	
	\code{ramwas1scanBams} executes \code{pipelineProcessBam} for each bam file.
	The bam processing is done faster using multiple CPU cores if \code{cputhreads} parameter is set to be greater than one.
	The function creates a local cluster with \code{param$cputhreads} nodes and runs the processing of BAMs on the cluster nodes.
}
\value{
	\code{pipelineProcessBam} returns a text string indicating the success (\code{"OK."}) or containing the error message.
	
	\code{ramwas1scanBams} returns a vector of success/error messages, one for each bam file.
}
\references{
	.
}
\seealso{
	See the vignettes.
}
\note{
	.
}
\author{
	Andrey Shabalin \email{ashabalin@vcu.edu}
}
\examples{
0
}
\keyword{ bam }
\keyword{ Rbam }
