---
bibliography: vignettes/ace.bib
---

# ACER

<!-- badges: start -->
<!-- badges: end -->

The goal of ACER is to test for differential essentiality
between sets of samples from a CRISPR knockout screen [@Hutton2020].  See
`vignettes(package=ACER)` for detailed examples.

## Installation

You can install the released version of ACER from GitHub with:

``` r
install_github("ACER")
```

## Basic Analysis

To determine differential essentiality, run the following commands (shown with example data):

``` r
library(ACER)
newDataObj <- DataObj$new(masterFiles = system.file('extdata','masterLibrarySeq.tsv', package='ACER'),
                          countFile = system.file('extdata','countData.tsv', package='ACER'),
                          negCtrlFile = system.file('extdata','negCtrlGenes.txt', package='ACER'), 
                          sampleInfoFile=system.file('extdata','sampleAnnotations.txt', package='ACER'),
                          hasInitSeq = T)
newModelObj <- ModelObj$new(user_DataObj = newDataObj,
                            use_neg_ctrl=T,
                            test_samples='test',
                            use_master_library = T)
newResultsObj <- optimizeModelParameters(user_DataObj = newDataObj,
                                         user_ModelObj = newModelObj)
writeResObj(newResultsObj)
```
# References
