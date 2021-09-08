---
---

# ACER

<!-- badges: start -->
[![Codecov test coverage](https://codecov.io/gh/CshlSiepelLab/ACE/branch/master/graph/badge.svg)](https://codecov.io/gh/CshlSiepelLab/ACE?branch=master)
[![Build Status](https://travis-ci.com/CshlSiepelLab/ACE.svg?token=ULtSN5KyvxhgFXqRqPas&branch=main)](https://travis-ci.com/CshlSiepelLab/ACE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5493071.svg)](https://doi.org/10.5281/zenodo.5493071)


<!-- badges: end -->

The goal of ACER is to test for differential essentiality
between sets of samples from a CRISPR knockout screen [[1]](#1).  See
`vignettes(package=ACER)` for detailed discussion.

## Installation

You can install the released version of ACER from GitHub with:

``` r
install_github("ACER")
```

## Basic Analysis

To determine differential essentiality, run the following commands (shown with example data):

``` r
library(ACER)
newDataObj <- DataObj$new(masterFiles = system.file('extdata','masterLibraryCounts.csv', package='ACER'),
                          countFile = system.file('extdata','countData.csv', package='ACER'),
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
<a id="1">[1]</a>
Hutton, E. R., Vakoc, C. R., and Siepel, A. (2020).
ACE: A Probabilistic Model for Characterizing Gene-Level Essentiality in CRISPR Screens.
bioRxiv, https://www.biorxiv.org/content/10.1101/868919v2
