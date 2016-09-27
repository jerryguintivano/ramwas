## ----generateData--------------------------------------------------------
{
	suppressPackageStartupMessages(library(ramwas))
	dir = paste0(tempdir(), '/simulated_project');
	ramwas0createArtificialData(dir, nsamples = 20, nreads = 1e6, verbose=FALSE)
	cat("Artificial project files created in:", "\n", dir, "\n")
	message("The generated files and directories are:");
	message(paste(list.files(dir), collate="\n"))
}

## ----parameters----------------------------------------------------------
param = list(
    dirproject = dir,
    dirbam = paste0(dir,"/bams"),
    scoretag = "mapq",
    minscore = 4,
    filebamlist = "bam_list.txt",
    filebam2sample =  "bam_list.txt",
    filecovariates = "covariates.txt",
    modeloutcome = "casecontrol",
    filecpgset = paste0(dir,"/Single_chromosome.rds"),
    minfragmentsize = 75,
    maxfragmentsize = 250,
    cputhreads = 2,
    minavgcpgcoverage = 0.3,
    minnonzerosamples = 0.3
)

## ----scan bams-----------------------------------------------------------
{
	suppressMessages({
	ramwas1scanBams(param)
	})
	message("The generated files and directories are:");
	message(paste(list.files(dir), collate="\n"))
}

## ----collectQC-----------------------------------------------------------
{
	suppressMessages({
	ramwas2collectqc(param)
	})
	message("The generated files and directories in QC directory are:");
	message(paste(list.files(paste0(dir,"/qc")), collate="\n"))
}

## ----normCoverage--------------------------------------------------------
{
	suppressMessages({
	ramwas3NormalizedCoverage(param)
	})
	message("The generated files and directories in coverage directory are:");
	message(paste(list.files(paste0(dir,"/coverage_norm_20")), collate="\n"))
}

## ----pca-----------------------------------------------------------------
{
	suppressMessages({
	ramwas4PCA(param)
	})
	message("The generated files and directories in PCA directory are:")
	message(paste(list.files(paste0(dir,"/coverage_norm_20/PCA_00_cvrts")), collate="\n"))
	message("Correlation of top PCs with covariates are:");
	head( read.table(paste0(dir,"/coverage_norm_20/PCA_00_cvrts/PC_vs_covariates_corr.txt"), sep = "\t", header = TRUE))
	message("P-values for the associuation of top PCs with covariates are:");
	head( read.table(paste0(dir,"/coverage_norm_20/PCA_00_cvrts/PC_vs_covariates_pvalue.txt"), sep = "\t", header = TRUE))
}

## ----mwas1---------------------------------------------------------------
{
	suppressMessages({
	param$modeloutcome = "casecontrol"
	ramwas5MWAS(param)
	})
	message("The generated files and directories in MWAS directories are:");
	message(paste(list.files(paste0(dir,"/coverage_norm_20/PCA_00_cvrts/Testing_casecontrol_0_PCs")), collate="\n"))
	message("The top MWAS findings are:");
	head( read.table(paste0(dir,"/coverage_norm_20/PCA_00_cvrts/Testing_casecontrol_0_PCs/Top_tests.txt"), sep = "\t", header = TRUE))
	qqPlotFast(fm.load(paste0(dir,"/coverage_norm_20/PCA_00_cvrts/Testing_casecontrol_0_PCs/Stats_and_pvalues"))[,3])
	title("Case-control status vs. methylation coverage")
}

## ----mwas2---------------------------------------------------------------
{
	suppressMessages({
	param$modeloutcome = "age"
	ramwas5MWAS(param)
	})
	message("The generated files and directories in MWAS directories are:");
	message(paste(list.files(paste0(dir,"/coverage_norm_20/PCA_00_cvrts/Testing_age_0_PCs")), collate="\n"))
	message("The top MWAS findings are:");
	head( read.table(paste0(dir,"/coverage_norm_20/PCA_00_cvrts/Testing_age_0_PCs/Top_tests.txt"), sep = "\t", header = TRUE))
	qqPlotFast(fm.load(paste0(dir,"/coverage_norm_20/PCA_00_cvrts/Testing_age_0_PCs/Stats_and_pvalues"))[,3])
	title("Age vs. methylation coverage")
}

## ----cv------------------------------------------------------------------
{
	suppressMessages({
	suppressWarnings({
	param$modeloutcomes = "age"
	ramwas6crossValidation(param)
	param$mmncpgs = 10
	ramwas7multiMarker(param)
	param$mmncpgs = 100
	ramwas7multiMarker(param)
	})
	})
	message("The generated files and directories in Multi-marker directory are:");
	message(paste(list.files(paste0(dir,"/coverage_norm_20/PCA_00_cvrts/Testing_age_0_PCs/CV_10_folds")), collate="\n"))
}

## ----clean---------------------------------------------------------------
unlink(dir, recursive=TRUE)

