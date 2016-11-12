\name{ramwas0createArtificialData}
\alias{ramwas0createArtificialData}
\title{
	Create Artificial Data Set
}
\description{
	Creates a set of artificial BAM files and supplementary files
	which can be used to test run the pipeline.
	The BAMs contain reads aligned only to one human chromosome, 
	with methylation effects embedded for simulated age and case-control status.
}
\usage{
ramwas0createArtificialData(dir, 
                            nsamples = 20, 
                            nreads = 1e6, 
                            ncpgs = 500e3, 
                            randseed = 18090212,
                            verbose = TRUE)
}
\arguments{
  \item{dir}{
    Directory for generated RaMWAS project files and BAMs.
  }
  \item{nsamples}{
    Number of sample/BAMs to create.
  }
  \item{nreads}{
    Number of reads in each BAM file.
  }
  \item{ncpgs}{
    Number of CpGs in the generated genome (with a single chromosome).
  }
  \item{randseed}{
    Random number generator seed for consistency of the output.
  }
  \item{verbose}{
    If TRUE, produce a message for each BAM file.
  }
}
\details{
  The function generates a number of files within \code{dir} directory.
  \enumerate{
     \item \code{bam_list.txt} - list of created BAM files. 
     To be used in \code{filebamlist} and 
     \code{filebam2sample} parameters in the pipeline.
     \item \code{covariates.txt} - table with age and sex status covariates. 
     For use in \code{filecovariates} parameter in the pipeline.
     \item \code{Single_chromosome.rds} - CpG location file 
     with the selected chromosome only.
     \item \code{bams} - directory with all the BAM files.
  }
}
\author{
	Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\seealso{
    For detailed description of all available parameters run:\cr 
    \code{vignette("RW5_param")}.
    
    For full description of the QC measures and the ploting options run\cr 
    \code{vignette("RW3_BAM_QCs")}
}
\examples{
### Location for the artificial project
dir = paste0(tempdir(), '/simulated_project');

ramwas0createArtificialData(dir, 
                            nsamples = 4, 
                            nreads = 10e3,
                            ncpgs = 1e3)
                            
message("Artificial project files created in:")
message(dir)
message("The generated files are:");
message(paste(list.files(dir, recursive=TRUE), collate="\n"))

### Clean up
unlink(dir, recursive=TRUE)
}
