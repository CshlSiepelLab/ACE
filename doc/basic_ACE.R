## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(ACER)

## ----sample inputs------------------------------------------------------------
head(fread(system.file('extdata','countData.csv', package='ACER')))

## ----basic--------------------------------------------------------------------
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

