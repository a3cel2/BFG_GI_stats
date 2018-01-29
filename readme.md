# Genetic interaction definitions for BFG-GI

Statistical analysis of genetic interactions mapped using Barcode Fusion Genetics

## Getting Started

To obtain the genetic interaction scores (GIS) for each gene pair 'xy' from the raw number of sequencing reads (C_xy) you just need to run 'master.R' . The only data used is table_s1.tsv (C_xy) and Cxy_Table_WithReplicates.txt in the 'data' folder. All newly made functions for master.R are in the packages folder. See below for other 'Dependencies'.  

We hexecute 'master.R' using the following steps:  
a) get into the ~/.../BFG_GI_stats-master/scripts directory  
b) Start an R session  
c) type source('master.R')

Calling master.R directly from R will work as well.

```
source('master.R')
```

To run the analysis on one of two counts instead of both replicates change this line in master.R

```
to_analyze <- 'sum'
```

to
```
to_analyze <- 'R1'

```
or
```
to_analyze <- 'R2'
```

## Folder Contents
**scripts/** - contains the main script (master.R), as well as development prototypes in the scratch folder

**data/** - contains:
1)  the input raw C_xy data (table_s1_input.tsv)
2)  output tables with GIS, Z_GIS and FDR at barcode-pair level (~/data/output/table_s1.tsv)
3)  output tables with GIS, Z_GIS and FDR at gene-pair level (~/data/output/table_s1_genewise.tsv)
4)  output tables with Delta_Z, Delta_GI and Delta_Z_FDR for each gene-pair for each condition-pair (~/data/output/table_s2_significant.tsv)
5)  Subset of (4) with significant pairs based on  Delta_Z_FDR (~/data/output/table_s2_significant.tsv)
Note - more outfiles are provided in 'results' folder

**doc/** - contains an explanation of the new genetic interaction definitions and the significance calling algorithm.  The latest documentation is included in the publication

**packages/** - contains all scripts that master.R depends on in the BfgGiStats subfolder, in an R package format

**results/** - all resulting plots from master.R


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

Cairo must also be properly configured to write PDF files
