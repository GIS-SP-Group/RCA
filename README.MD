# RCA version 2.0
RCA, (Reference Component Analysis), is a computational approach for robust cell type annotation of single cell RNA sequencing data (scRNAseq). 
It is developed by the Prabhakar lab at the Genome Institute of Singapore (GIS). 
The original version of RCa was published in Nature Genetics (doi: 10.1038/ng.3818, Li et al., 2017).

## Release notes
Version 2.0 
Release date: September 3, 2019

## Functionality of RCA
Features: 
* clustering analysis of scRNAseq data from Human samples
* two modes :
  * "GlobalPanel": default option for clusterig general single cell data sets that include a wide spectrum of cell types.          
  * "ColonEpitheliumPanel": suitable for analyzing human colon/intestine derived samples.


## A beginner's guide to RCA
This guide will walk you through installing RCA and will showcase a exemplary analysis of publicly available scRNA-seq data from *10X Genomics*.

### Install the RCA R-package
Before you try to install RCA make sure that your R-version is at least R 3.5.0.
You can directly install RCA from github using the commands:

```{r}
library(remotes)
install_github("prabhakarlab/RCAv2")
```

The current release of RCA requires the following packages to be available on your system:
* remotes
* Matrix
* qlcMatrix
* WGCNA
* fastcluster
* BiocManager
* ggplot2
* gridExtra
* dplyr
* ComplexHeatmap
* circlize
* umap
* ggpubr

All missing CRAN-packages can be automatically installed during the RCA installation.
Note that the *HiClimR* package is optional, it allows to speed up the computationally expensive steps.

### Load the package
After installation, load the package with the command
```{r}
library("RCAv2")
```

### Generate a RCA object from scRNA-seq data



### Compute a projection to a (custom) reference data set


## FAQ
What is the difference between the original RCA version (Li et al., Nat Genet, 2017) and RCA version 2?
*Answer here*

The clustering is very slow (or it doesn't work at all), what can I do?
*Answer here*

The cell-type I am interested in are not part of the provided panels, can I generate and use my own reference data set in RCA?
*Answer here*

Why is RCA currently available only in R?
*Answer here*

I do not have 10X data, can I just use a count matrix as input for RCA?
*Answer here*

How to cite RCA?
*Answer here*
