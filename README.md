# RaMWAS: Fast methylome-wide association study pipeline for enrichment platforms

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
source("https://bioconductor.org/biocLite.R")
biocLite("ramwas")
```

### Update to GitHub Version

To update RaMWAS to the development version from GitHub, run

```
devtools::install_github("andreyshabalin/ramwas")
```

If `devtools` package is missing, it can be installed with

```
install.packages("devtools")
```

### Run Rcmd check

To rerun package tests, suggested packages can be installed with

```
source("https://bioconductor.org/biocLite.R")
biocLite(c("BSgenome.Ecoli.NCBI.20080805","knitr","BiocCheck", 
            "rmarkdown", "pander", "BiocStyle", "devtools"))
```
