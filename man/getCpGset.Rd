\name{getCpGset}
\alias{getCpGset}
\title{
	Construct CpG set for a Reference Genome
}
\description{
Finds all CpGs in a reference genome.
}
\usage{
getCpGset(genome)
}
\arguments{
  \item{genome}{
    A \code{\linkS4class{BSgenome}} object or a character vector with genome sequences.
  }
}
\details{
  The function searches for all \code{CG} pairs in the genome.
}
\value{
	Returns a list with CpG coordinates for each genome sequence.
}
\author{
	Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\examples{
library(BSgenome.Ecoli.NCBI.20080805)
cpgset = getCpGset(BSgenome.Ecoli.NCBI.20080805);

print('First 10 CpGs in NC_008253:')
print(cpgset$NC_008253[1:10])


genome = list(
  chr1 = 'AGCGTTTTCATTCTGACTGCAACGGGC',
  chr2 = 'AAAAAACGCCTTAGTAAGTGATTTTCG');
cpgset = getCpGset(genome);

print('CpG coordinates in the toy genome:')
print(cpgset)
}
