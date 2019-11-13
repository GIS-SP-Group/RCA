# RCA version 2.0
RCA, short for Reference Component Analysis, is a R package for robust cell type annotation of single cell RNA sequencing data (scRNAseq). It is developed by Prabhakar Group at the Genome Institute of Singapore (GIS).

## Versions and Features

Version 2.0 
Release date: September 3, 2019

Features: 
* clustering analysis of scRNAseq data from Human samples
* two modes :
  * "GlobalPanel": default option for clusterig general single cell data sets that include a wide spectrum of cell types.          
  * "ColonEpitheliumPanel": suitable for analyzing human colon/intestine derived samples.


## A beginner's guide to RCA
### Install the RCA R-package
You can directly install RCA from github using the commands:

```{r}
library(remotes)
install_github("prabhakarlab/RCAv2")
```
