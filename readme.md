# Genetic interaction definitions for BFG-GI

Statistical analysis of genetic interactions mapped using Barcode Fusion Genetics (BFG-GI)

## Getting Started

To analyze genetic interactions from the raw number of sequencing reads (C_xy) run 'master.R' . The input data used is table_s2.txt. All functions called by master.R are in the packages folder. See below for other 'Dependencies'.  This script performs several tasks - including GIS calculation, Z_GIS calculation, FDR calculation, genetic interaction calls, and differential genetic interaction calls. It reproduces several plots used in the manuscript

We execute 'master.R' using the following steps:  
a) get into the ~/.../BFG_GI_stats-master/scripts directory  
b) Start an R session  
c) type source('master.R')

Calling master.R directly from R will work as well.  To properly find the input data, the working directory must be the same as the master.R script.

```
source('master.R')
```

## Folder Contents
**scripts/** - contains the main script (master.R)

**data/** - contains:
1)  the input data (table_s2.txt)
2)  tables resulting from analysis in the output/ folder - tables s3,s4,s5,s6 in the publication

**doc/** - contains an explanation of the new genetic interaction definitions and the significance calling algorithm.  The latest documentation is included in the publication

**packages/** - contains all scripts that master.R depends on in the BfgGiStats subfolder, in an R package format

**results/** - all plots resulting from from master.R


## Authors

**Albi Celaj** (https://github.com/a3cel2)  
**Javier Diaz (testing and documentation)** (https://github.com/jdime)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## Dependencies
The following R packages must be installed prior to running master.R:
1) dplyr
2) devtools
3) roxygen2
4) ROCR
5) ExtDist
6) qvalue (Bioconductor)

Cairo must be installed for accurate reproduction of PDF files.  Otherwise, the base pdf funcion will be used in master.R
