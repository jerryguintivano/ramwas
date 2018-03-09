\name{injectSNPsMAF}
\alias{injectSNPsMAF}
\title{
    Inject SNPs from VCF Count File into a DNA Sequence
}
\description{
    Injects SNPs from a VCF count file into a DNA sequence.
}
\usage{
injectSNPsMAF(gensequence, frqcount, MAF = 0.01)
}
\arguments{
    \item{gensequence}{
        A string or \code{\linkS4class{DNAString}} of the DNA sequence.
    }
    \item{frqcount}{
        File name of the allele count file produced by
    \code{vcftools} with \code{--counts} parameter.
        Alternatively, the file content can be provided as a
        character vector (see \code{\link[base]{readLines}}).
    }
    \item{MAF}{
        SNPs with minor allele frequency at or above \code{MAF} are injected.
    }
}
\value{
    Returns a string with the genome sequence with SNPs injected.
}
\author{
    Andrey A Shabalin \email{andrey.shabalin@gmail.com}
}
\seealso{
    See \code{\link[BSgenome]{injectSNPs}} for the standard analog
    function without MAF filtering.
}
\examples{
gensequence1 = "AAAACAAAA"
frqcount = c(
    "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:COUNT}",
    "1\t6\t2\t1000\tA:400\tG:600",
    "1\t7\t2\t1000\tA:800\tC:200",
    "1\t9\t2\t1000\tA:900\tG:100")
MAF = 0.01

gensequence2 = injectSNPsMAF(gensequence1, frqcount, MAF)

### No CpGs without SNPs
show(gensequence1)
getCpGsetCG(gensequence1)

### SNPs create 1 CpG
show(gensequence2)
getCpGsetALL(gensequence2)
}
