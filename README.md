# RaMWAS: Fast methylome-wide association study pipeline for enrichment platforms

RaMWAS provides a complete toolset for 
    methylome-wide association studies (MWAS).
    It is specifically designed for data from 
    enrichment based methylation assays,
    but can be applied to other methylomic data as well.
    The analysis pipeline includes seven steps:  
    (1) scanning aligned reads from BAM files,  
    (2) calculation of quality control measures,  
    (3) creation of methylation score (coverage) matrix,  
    (4) principal component analysis for capturing batch effects
    and detection of outliers,  
    (5) association analysis with respect to phenotypes of interest
    while correcting for top PCs and known covariates,  
    (6) annotation of significant findings, and  
    (7) multi-marker analysis (methylation risk score) using elastic net.  
    Additionally, RaMWAS include tools for joint analysis of methlyation and
    genotype data.
    
------------
INSTALLATION
------------

## Install via Bioconductor

To install RaMWAS via Bioconductor run:

```
source("http://bioconductor.org/biocLite.R")
biocLite("ramwas")
```

## Update via GitHub

To update RaMWAS to the development version from GitHub run

```
devtools::install_github("andreyshabalin/ramwas")
```

If the `devtools` package is not installed, run this line prior:

```
install.packages("devtools")
```
