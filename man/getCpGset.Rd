\name{getCpGset}
\alias{getCpGsetCG}
\alias{getCpGsetALL}
\title{
    Construct CpG set for a Reference Genome
}
\description{
    Finds all CpGs in a reference genome.
}
\usage{
getCpGsetCG(genome)
getCpGsetALL(genome)
}
\arguments{
    \item{genome}{
        A \code{\link[BSgenome]{BSgenome}} object or a
        character vector with genome sequences.
    }
}
\details{
    The \code{getCpGsetCG} function searches for all
    \code{CG} pairs in the genome. \cr
    The \code{getCpGsetALL} function also works for genomes with injected SNPs.
}
\value{
    Returns a list with CpG coordinates for each genome sequence.
}
\author{
    Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\examples{
### Using a BSGenome input

library(BSgenome.Ecoli.NCBI.20080805)
cpgset = getCpGsetCG(BSgenome.Ecoli.NCBI.20080805)

print("First 10 CpGs in NC_008253:")
print(cpgset$NC_008253[1:10])


### Using a character vector input

genome = list(
    chr1 = "AGCGTTTTCATTCTGACTGCAACGGGCYR",
    chr2 = "AAAAAACGCCTTAGTAAGTGATTTTCGYR")
cpgset1 = getCpGsetCG(genome)
cpgset2 = getCpGsetALL(genome)

print("Pure CG coordinates in the toy genome:")
print(cpgset1)

print("CG coordinates in the toy genome possible with SNPs:")
print(cpgset2)
}
