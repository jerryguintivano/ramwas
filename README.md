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

## Install from GitHub

### Prerequisites

To install RaMWAS several R packages must be installed first

```
install.packages(c("knitr","rmarkdown","KernSmooth","filematrix","digest","glmnet","devtools","pander"))
source("https://bioconductor.org/biocLite.R")
biocLite(c("BiocInstaller","BiocStyle","GenomicAlignments","Rsamtools","biomaRt"))
```

### Installation

To install RaMWAS directly from GitHub run

```
devtools::install_github("andreyshabalin/ramwas")
```

### Prerequisites for building vignettes

The following packages are used in vignettes only

```
source("https://bioconductor.org/biocLite.R")
biocLite(c("BSgenome.Hsapiens.UCSC.hg19","SNPlocs.Hsapiens.dbSNP144.GRCh37","BSgenome.Ecoli.NCBI.20080805"))
```
