\name{insilicoFASTQ}
\alias{insilicoFASTQ}
\title{
    Construct FASTQ File for In-silico Alignment Experiment
}
\description{
    Creates a FASTQ file with all fragments of \code{fraglength} bp long.
}
\usage{
insilicoFASTQ(con, gensequence, fraglength)
}
\arguments{
    \item{con}{
        A \code{\link[base:connections]{connection}} object or a
        character string naming the output file. \cr
        If the name ends with ".gz", a compressed file is created. \cr
        An empty string can be used to output to the console.
    }
    \item{gensequence}{
        A string or \code{\linkS4class{DNAString}} of the DNA sequence.
    }
    \item{fraglength}{
        Fragment length.
    }
}
\details{
    The function a FASTQ file with all fragments of \code{fraglength} bp
    long from the forward strand of the DNA sequence.
}
\value{
    Returns a list with CpG coordinates for each genome sequence.
}
\author{
    Andrey A Shabalin \email{andrey.shabalin@gmail.com}
}
\examples{
## There are four 4 bp fragments in a 7 basepair sequence:
insilicoFASTQ(con="", gensequence = "ABCDEFG", fraglength=4)
}
