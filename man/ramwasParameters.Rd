\name{ramwasParameters}
\alias{ramwasParameters}
\title{
    Function for Convenient Filling of the RaMWAS Parameter List.
}
\description{
    RaMWAS parameter vector which is used by major functions of the pipeline
    is a regular R list and setting it does not require a special function.
    However, using this function makes it much simpler in RStudio
    as the names and role of every parameter is showed in the RStudio IDE.
}
\usage{
ramwasParameters(
    dirproject,
    dirfilter,
    dirrbam,
    dirrqc,
    dirqc,
    dircoveragenorm,
    dirtemp,
    dirpca,
    dirmwas,
    dircv,
    dirbam,
    filebamlist,
    bamnames,
    filebam2sample,
    bam2sample,
    filecpgset,
    filenoncpgset,
    filecovariates,
    covariates,
    cputhreads,
    diskthreads,
    usefilelock,
    scoretag,
    minscore,
    maxrepeats,
    minavgcpgcoverage,
    minnonzerosamples,
    buffersize,
    doublesize,
    modelcovariates,
    modeloutcome,
    modelPCs,
    modelhasconstant,
    qqplottitle,
    toppvthreshold,
    mmncpgs,
    mmalpha,
    cvnfolds,
    bihost,
    bimart,
    bidataset,
    biattributes,
    bifilters,
    biflank,
    fileSNPs,
    dirSNPs,
    ...)
}
\arguments{
    \item{dirproject}{
        The project directory. Default is currect directory. \cr
        Files specified by "file*" parameters are looked for here,
        unless they have full path specified.
    }
    \item{dirfilter}{
        By default, the same as "dirproject". \cr
        All files created by RaMWAS are created within this directory. \cr
        If the user wants to test different read filtering rules,
        they can dirfilter to TRUE.
        This will set it to something like "Filter_MAPQ_4",
        there "MAPQ" is the BAM field used for filtering and
        "4" is the thredhold.
    }
    \item{dirrbam}{
        Directory where RaMWAS saves RaMWAS raw data files
        (read start locations) after scanning BAMs. \cr
        It is "rds_rbam" by default and located in "dirfilter".
    }
    \item{dirrqc}{
        Directory where RaMWAS saves QC files in R format after scanning BAMs.
        \cr
        It is "rds_qc" by default and located in "dirfilter".
    }
    \item{dirqc}{
        Directory where RaMWAS saves QC plots and text files (BAM QC info)
        after scanning BAMs. \cr
        It is "qc" by default and located in "dirfilter".
    }
    \item{dircoveragenorm}{
        Directory where RaMWAS saves coverage matrix at Step 3 of the pipeline.
        \cr
        It is "coverage_norm_123" by default (123 is the number of samples)
        and located in "dirfilter".
    }
    \item{dirtemp}{
        Directory where RaMWAS stores temporary files
        during construction of coverage matrix at Step 3 of the pipeline. \cr
        It is "temp" by default and located in "dircoveragenorm". \cr
        For better performance it can be set to a location
        on a different hard drive than "dircoveragenorm".
    }
    \item{dirpca}{
        Directory where RaMWAS saves results of PCA analysis at Step 4
        of the pipeline. \cr
        It is "PCA_12_cvrts_0b0a0c" by default and located in "dircoveragenorm",
        where 12 is the number of covariates regressed out and
        "0b0a0c" is a unique code to differenciate
        different sets of 12 covariates.
    }
    \item{dirmwas}{
        Directory where RaMWAS saves results of MWAS analysis at Step 5
        of the pipeline. \cr
        It is "Testing_age_7_PCs" by default and located in "dirpca", where
        "age" is the phenotype being tested and
        "7" is number of top PCs included in the model.
    }
    \item{dircv}{
        Directory where RaMWAS saves results of Methylation Risk Score analysis
        at Step 7 of the pipeline. \cr
        It is "CV_10_folds" by default and located in "dirmwas", where
        10 is number of folds in N-fold cross validation.
    }
    \item{dirbam}{
        Location of BAM files. \cr
        If not absolute, it is considered to be relative to "dirproject".
    }
    \item{filebamlist}{
        If defined, must point to a text file with one BAM file name per line.
        \cr
        BAM file names may include path, relative to "dirbam" or absolute. \cr
    }
    \item{bamnames}{
        A character vector with BAM file names. \cr
        Not required if "filebamlist" is specified. \cr
        BAM file names may include path, relative to "dirbam" or absolute.
    }
    \item{filebam2sample}{
        Allowes multiple BAMs contain information about common sample. \cr
        Must point to a file with lines like "sample1=bam1,bam2,bam3".
    }
    \item{bam2sample}{
        Allowes multiple BAMs contain information about common sample. \cr
        Not required if "filebam2sample" is specified. \cr
        Must be a list like \code{list(sample1 = c("bam1","bam2","bam3"),
        sample2 = "bam2")}
    }
    \item{filecpgset}{
        Name of the file storing a set of CpGs.
    }
    \item{filenoncpgset}{
        If defined, must point to a file storing vetted
        locations away from any CpGs.
    }
    \item{filecovariates}{
        Name of the file containing phenotype and covariates
        for the available samples. \cr
        If the file has extension ".csv",
        it is assumed to be comma separated,
        otherwise - tab separated.
    }
    \item{covariates}{
        Data frame with phenotype and covariates
        for the available samples. \cr
        Not required if "filecovariates" is specified.
    }
    \item{cputhreads}{
        Maximum number of CPU intensive tasks running in parallel. \cr
        Set to the number of CPU cores by default.
    }
    \item{diskthreads}{
        Maximum number of disk intensive tasks running in parallel. \cr
        Set to 2 by default.
    }
    \item{usefilelock}{
        If TRUE, parallel jobs are prevented from simultaneous access
        to file matrices. \cr
        Can improve performance on some systems.
    }
    \item{scoretag}{
        Reads from BAM files are filtered by this tag. \cr
        The "minscore" parameter defines the minimum admissible score.
    }
    \item{minscore}{
        Reads from BAM files with score "scoretag" below this are excluded.
    }
    \item{maxrepeats}{
        Duplicate reads (reads with the same start position and direction)
        in excess of this limit are removed.
    }
    \item{minavgcpgcoverage}{
        CpGs with average coverage below this threshold are removed.
    }
    \item{minnonzerosamples}{
        CpGs with fraction of samples with non-zero coverage below this
        threshold are removed.
    }
    \item{buffersize}{
        Coverage matrix transposition is performed using buffers of this size.
        \cr
        Larger "buffersize" improves speed of Step 3 of the pipeline, but
        requires more memory. \cr
        Default is 1e9, i.e. 1 GB.
    }
    \item{doublesize}{
        The coverage matrix is stored with this number of bytes per value. \cr
        Set to 8 for full (double) precision. \cr
        Set to 4 to use single precision and
        create 50\% smaller coverage filematrix.
    }
    \item{modelcovariates}{
        Names of covariates included in PCA and MWAS.
    }
    \item{modeloutcome}{
        Name of the outcome variable for MWAS.
    }
    \item{modelPCs}{
        Number of principal components accounted for in MWAS.
    }
    \item{modelhasconstant}{
        By default, the tested linear model includes a constant. \cr
        To exclude it, set "modelhasconstant" parameter to \code{FALSE}.
    }
    \item{qqplottitle}{
        The title of the QQ-plot produced by MWAS (step 4 of the pipeline).
    }
    \item{toppvthreshold}{
        Determines the number of top MWAS results saved in text file. \cr
        If it is 1 or smaller, it defines the p-value threshold. \cr
        If larger than 1, it defines the exact number of top results.
    }
    \item{mmncpgs}{
        Parameter for multi-marker elastic net cross validation (MRS). \cr
        Defines the number of top CpGs on which to train the elastic net. \cr
        Can be set of a vector of multiple values, each is tested separately.
    }
    \item{mmalpha}{
        Parameter for multi-marker elastic net cross validation (MRS). \cr
        Elastic net mixing parameter alpha. \cr
        Set to 0 by default.
    }
    \item{cvnfolds}{
        Parameter for multi-marker elastic net cross validation (MRS). \cr
        The number of folds in the N-fold cross validation.
    }
    \item{bihost}{
        Parameter for BiomaRt annotation (Step 6 of the pipeline). \cr
        BioMart host site. \cr
        Set to \code{"grch37.ensembl.org"} by default.
    }
    \item{bimart}{
        Parameter for BiomaRt annotation (Step 6 of the pipeline). \cr
        BioMart database name, see \code{\link[biomaRt]{listMarts}}. \cr
        Set to \code{"ENSEMBL_MART_ENSEMBL"} by default.
    }
    \item{bidataset}{
        Parameter for BiomaRt annotation (Step 6 of the pipeline). \cr
        BioMart data set, see \code{\link[biomaRt]{listDatasets}}. \cr
        Set to \code{"hsapiens_gene_ensembl"} by default.
    }
    \item{biattributes}{
        Parameter for BiomaRt annotation (Step 6 of the pipeline). \cr
        BioMart attributes of interest,
        see \code{\link[biomaRt]{listAttributes}}. \cr
        Set to \code{c("hgnc_symbol","entrezgene","strand")} by default.
    }
    \item{bifilters}{
        Parameter for BiomaRt annotation (Step 6 of the pipeline). \cr
        BioMart filters (if any), see \code{\link[biomaRt]{listFilters}}. \cr
        Set to \code{list(with_hgnc_transcript_name=TRUE)} by default ignore
        genes without names.
    }
    \item{biflank}{
        Parameter for BiomaRt annotation (Step 6 of the pipeline). \cr
        Allowed distance between CpGs and genes
        or other annotation track elements. \cr
        Set to 0 by default, requiring direct overlap.
    }
    \item{fileSNPs}{
        Name of the filematrix with genotype (SNP) data. \cr
        The filematrix dimensions must match the coverage matrix.
    }
    \item{dirSNPs}{
        Directory where RaMWAS saves the results of joint
        methylation-genotype analysis.
    }
    \item{\dots}{
        Any other named parameters can be added here.
    }
}
\details{
    The function simply collects all the parameters in a list. \cr
    The main benefit of the function is that the user does not
    need to memorize the names of RaMWAS parameters.

    Here is how it helps in RStudio: \figure{hint.png}
}
\value{
    List with provided parameters.
}
\author{
    Andrey A Shabalin \email{andrey.shabalin@gmail.com}
}
\seealso{
    See vignettes: \code{browseVignettes("ramwas")}.
}
\examples{
ramwasParameters(dirproject = ".", cputhreads = 4)
}
