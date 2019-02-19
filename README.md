# RaMWAS: Fast Methylome-Wide Association Study Pipeline for Enrichment Platforms

RaMWAS provides a complete toolset for 
methylome-wide association studies (MWAS).
It is specifically designed for data from 
enrichment based methylation assays,
but can be applied to other methylomic data as well.
The analysis pipeline includes seven steps:

* scanning aligned reads from BAM files,  
* calculation of quality control measures,  
* creation of methylation score (coverage) matrix,  
* principal component analysis for capturing batch effects
and detection of outliers,  
* association analysis with respect to phenotypes of interest
while correcting for top PCs and known covariates,  
* annotation of significant findings, and  
* multi-marker analysis (methylation risk score) using elastic net.  

Additionally, RaMWAS include tools for joint analysis of methlyation and
genotype data.

## Installation

### Install Bioconductor Version

To install
[Bioconductor version](https://bioconductor.org/packages/ramwas/)
of RaMWAS, run

```
if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ramwas")
```

### Update to GitHub Version

To update RaMWAS to the development version from GitHub, run

```
if (!requireNamespace("devtools", quietly = TRUE))
   install.packages("devtools")
devtools::install_github("andreyshabalin/ramwas")
```

### Run Rcmd Check (Install Suggested Packages)

To rerun package tests, suggested packages can be installed with

```
if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("BSgenome.Ecoli.NCBI.20080805", "knitr", "BiocCheck",
                       "rmarkdown", "pander", "BiocStyle", "devtools"))
```
