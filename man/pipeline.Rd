\name{pipeline}
\alias{ramwas1scanBams}
\alias{pipelineProcessBam}
\alias{ramwas2collectqc}
\alias{ramwas3NormalizedCoverage}
\alias{ramwas4PCA}
\alias{ramwas5MWAS}
\alias{ramwas6crossValidation}
\alias{ramwas7multiMarker}

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

ramwas3NormalizedCoverage(param)

ramwas4PCA(param)

ramwas5MWAS(param)

ramwas6crossValidation(param)

ramwas7multiMarker(param)
}
\arguments{
  \item{param}{
    List with RaMWAS parameters. The detailed description of all available parameters see at: \link[ramwas]{parameters}.
  }
  \item{bamname}{
    Name of the BAM file to process. Can be absolute or relative to \code{dirbam} parameter.
  }
}
\details{

	The RaMWAS pipeline consists of 7 major components.
	\enumerate{
		\item \bold{Scan BAM files, filtering reads, calculation of QC measures, saving the read start positions.} \cr 
		Function \code{ramwas1scanBams} scans all BAMs, while \code{pipelineProcessBam} can be used to process a single BAM.
		
		\itemize{
			\item BAM file names are taken from either \var{filebamlist} file or \var{bamnames} variable.
			\item BAMs can be located in different directories, but must have unique names.
			\item BAMs are looked for in \var{dirbam} directory.
			\item Reads are filtered by \var{scoretag} tag in the BAM file (typically \code{"MAPQ"} or \code{"AS"}). Reads with the score tag below \var{minscore} are excluded.
			\item Up to \var{cputhreads} \pkg{parallel} R jobs are used to process BAMs (if \var{cputhreads} > 1).
			\item The QC information plots are recorded in \var{dirqc} directory, which is \code{"qc"} by default.
			\item The QC information in R format is saved in \var{dirrqc} directory, which is \code{"rds_qc"} by default.
			\item The read start locations are saved in \var{dirrbam} directory, which is \code{"rds_rbam"} by default.
			\item If a BAM has already been scanned, it is not processed again.
			\item The QC (Quality Control) measures are calculated using a number of parameters, including \var{filecpgset}, \var{filenoncpgset}, \var{maxfragmentsize}. For more on QC measures see: \link[ramwas]{QCs}.
		}
		
		\item \bold{Collect QC measures from individual BAMs, summarizing by samples, recording summary.} \cr
		Function \code{ramwas2collectqc} collects QC information from individual scanned BAM files and created additional QC plots and files.
		\itemize{
			\item Multiple BAMs can contain information from a single sample. The bam to sample matching information is taken from \var{filebam2sample} file or directly from \var{bam2sample} character vector.
			\item All QC measures are summarized by BAM, by sample, and combined in one total metric across all bams. 
			\item The resulting plots and tables are recored in \var{dirqc} subdirectories \code{"summary_bams"}, \code{"summary_bams_in_bam2sample"}, \code{"summary_by_sample"}, and \code{"summary_total"}.
			\item The fragment size distribution is estimated using all BAMs in \var{bam2sample} (if defined) and saved in \code{"Fragment_size_distribution.txt"} file in project directory.
		}
		
		\item \bold{Calculate coverage matrix for all CpGs across all samples, filter CpGs by coverage, normalize samples to the same sample average.} \cr
		Function \code{ramwas3NormalizedCoverage} performs this part.
		\itemize{
			\item The coverage matrix is calculated for the samples listed in \var{bam2sample}. If \var{filecovariates} or \var{covariates} is defined, only the samples in the covariate file are used.
			\item The coverage matrix is created in \var{dircoveragenorm} directory. It is \code{"coverage_norm_1234"} by default, where 1234 is the number of samples.
			\item The coverage matrix is stored as a \code{\link[filematrix]{filematrix}} named \code{"Coverage"}. Filematrices can be accessed (by parts) as usual R matrices without the need to load the whole filematrix in memory.
			\item The coverage is calculated for CpGs in \var{filecpgset}.
			\item CpGs with less than \var{minavgcpgcoverage} average coverage are excluded (dafault is 0.3).
			\item CpGs with fewer than \var{minnonzerosamples} samples with nonzero coverage are excluded (default is 0.3, meaning a CpG is preserved if at least 30\% of samples must have non-zero coverage).
			\item The locations of CpG that passed the filter are saved in filematrix \code{"CpG_locations"}, which has two columns indicating chromosome number and position. The chomosome are numbered as in \var{filecpgset} (and chromosome names are also saved in \code{"CpG_chromosome_names.txt"}.
			\item Up to \var{cputhreads} \pkg{parallel} R jobs are used for calculation of raw coverage (This is a CPU intensive task).
			\item Up to \var{diskthreads} \pkg{parallel} R jobs are used for CpG filtering and sample normalization  (This is a disk heavy task).
			\item The temporary files are located in \var{dirtemp} directory, which is \code{"temp"} subdirectory of \var{dircoveragenorm} by default. \cr
				For better performance, set \var{dirtemp} to a location on another hard drive.
		}
		\item \bold{Principal Component Analysis} \cr
		Function \code{ramwas4PCA} runs PCA.
		\itemize{
			\item The principal components are calculated after removing the effect of selected covariates from the data. The covariates can be listed by names \var{modelcovariates} parameter. The covariates are taken from \var{filecovariates} or \var{covariates}.
			\item The PCA output is saved in \var{dirpca} directory, by default is named like \code{"PCA_05_cvrts_2d33b43"} and located in \var{dircoveragenorm}. The name includes the number of covariates excluded. To differenciate between different sets of, say, 5 covariates it also includes an 8 letter code, which is different for different sets of covariates (with a very high probability).
			\item Up to \var{diskthreads} \pkg{parallel} R jobs are used PCA.
			\item The sample covariance matrix is saved in \code{"covmat.rds"} file.
			\item The eigenvalue decomposition of the sample covariance matrix is saved in \code{"eigen.rds"} file.
			\item PC loadings are saved in text format in \code{"PC_loadings.txt"}.
			\item PCs are saved in text file \code{"PC_values.txt"}.
			\item Up to top 40 PCs are plotted in \code{"PC_plot_covariates_removed.pdf"}. The pdf also contains plots of PC loadings for the first few PCs.
			\item The correlations of PCs with covariates is recorded in \code{"PC_vs_covariates_corr.txt"}. The statistical significance of these correlations is recorded in \code{"PC_vs_covariates_pvalue.txt"}. \cr 
			P-values in \code{"PC_vs_covariates_pvalue.txt"} can be helpful in selection of covariates for MWAS (Methylome-wide association study). 
		}
		\item \bold{Methylome-wide association study (MWAS)} \cr
		Function \code{ramwas5MWAS} tests for association between normalized CpG coverage and a phenotype (covariate of interest) while accounting for a set of covariates and top principal components.
		\itemize{
			\item The phenotype of interest is named by \var{modeloutcome} parameter (taked from \var{covariates}).
			\item As in PCA, \var{modelcovariates} lists covariates to correct for in MWAS.
			\item The analysis can also include several top PCs as covariates. The number of PCs is set by \var{modelPCs}.
			\item The testing is performed using linear regression model.
			\item The MWAS output is saved in \var{dirmwas} directory, by default named like \code{"Testing_Age_2_PCs"} where "Age" with the phenotype of interest and 2 is the number of PCs included as covariates.
			\item Up to \var{diskthreads} \pkg{parallel} R jobs are used MWAS.
			\item The t-statistics, p-values and q-values for all CpGs are recorded in \code{"Stats_and_pvalues"} filematrix. The row of this filematrix match those of \code{"CpG_locations"}.
			\item The QQ-plot of all p-values, with an estimate of inflation parameter lambda is saved in \code{"QQ_plot.pdf"}. \cr
			The title of the qqplot can be set with the parameter \var{qqplottitle}.
		}
		\item \bold{Cross Validation with multi-marker elastic net model.} \cr
		RaMWAS pipeline can estimate predictive power of elastic net applied to top CpGs from MWAS. \cr
		To have an unbiased estimate of predictive power we use cross validation technique. \cr
		First, \code{ramwas6crossValidation} runs MWAS 10 times, each time excluding different 10\% of the data. \cr
		Then, \code{ramwas7multiMarker} trains elastic net on the top CpGs from each of 10 MWASes and uses it to predict the phenotype of interest on the remaing 10\%.
		
				Namely, the elastic net is trained on 90\% of the data and then used to predic the phenotype of interest on the remaing 10\%. This procedure is then applied 10 times to have a prediction for all samples.
		\itemize{
			\item By default, 10-fold cross validatidation is performed. The number of folds can be set by \var{cvnfolds} parameter.
			\item The split of samples into folds is random. For consistency across runs, the random seed can be set by \var{randseed}. It is set to 18090212 by default.
			\item The number of CpGs used to by elastic net is set by \var{mmncpgs} parameter. It is 1,000 by default.
			\item The elastic net parameter labmda can be set by \var{mmalpha} parameter. It is 0 by default.
		}

		
	}
}
\value{modelcovariates
	\code{pipelineProcessBam} returns "OK." if ran successfully, an error message otherwise.
	Other functions return \code{NULL}, all the output is saved in the respective directories.
}
%\references{
%}
\author{
	Andrey A Shabalin \email{ashabalin@vcu.edu}
}
\seealso{
For more on QC measures see: \link[ramwas]{QCs}. \cr
For more on pipeline parameters see: \link[ramwas]{parameters}.
}
\examples{
library(ramwas)
param = list(
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
ramwas1scanBams( param )

ramwas2collectqc( param )

ramwas3NormalizedCoverage(param)

ramwas4PCA(param)

ramwas5MWAS(param)

ramwas6crossValidation(param)

ramwas7multiMarker(param)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
