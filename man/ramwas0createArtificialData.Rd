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
        Number of samples/BAMs to create.
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

    The generated BAMs have 600 CpGs affected by \code{sex},
    namely fully methylated or not methylated at all, depending on \code{sex}.
    The methylation level of 1\% of all CpGs is affected by age.
    The methylation of those CpGs is
    equal to \code{age/100} or \code{1-age/100}.
    The age is generated randomly in the range from 20 to 80.
}
\value{
    The function creates multiple files but returns no value.
}
\author{
    Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\seealso{
    See vignettes: \code{browseVignettes("ramwas")}.
}
\examples{
### Location for the artificial project
dr = paste0(tempdir(), "/simulated_project")

ramwas0createArtificialData(
    dr,
    nsamples = 4,
    nreads = 10e3,
    ncpgs = 1e3)

# Artificial project files created in:
dr
# The generated files are:"
as.matrix(list.files(dr, recursive=TRUE))

### Clean up
unlink(paste0(dr,"/*"), recursive=TRUE)
}
