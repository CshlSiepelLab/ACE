context('Test Model Object')
library(ACER)

test_that('Master library replicate handling in getMasterLib.', {
  testMasterCountList <- list('masterA' = data.table('sgrna'=paste0('g', 1:1e4),
                                                     'rep1' = rep(0:9, 1e3)) ,
                              'masterB' = data.table('sgrna' = paste0('g', 1:1e4),
                                                     'repli1' = 1:1e4,
                                                     'repli2' = (1+2e4):3e4),
                              'masterC' = data.table('sgrna' = paste0('g', 1:1e3),
                                                     'rep1' = rep(0, 1e3)))
  print(head(testMasterCountList[[1]]))
  capture.output(testMasterFreq <- getMasterLib(countList = testMasterCountList),
                 file = 'getMasterLib_debugger.txt')
  noInfValuesInMasterFreq <- !testMasterFreq[, any(sapply(.SD, is.infinite)),
                                             .SDcols = -'sgrna']
  expect_true(noInfValuesInMasterFreq)
  noNaValuesInMasterFreq <- !testMasterFreq[, any(sapply(.SD, is.na)),
                                            .SDcols = -'sgrna']
  expect_true(noNaValuesInMasterFreq)
  expect_true(is.data.table(testMasterFreq))
})

test_that('ModelObj with blank guides has problem.', {
  if (!'testCountOne.txt' %in% dir()) {
    write.table(data.table('guide' = paste0('g', 1:3),
                           'gene' = 'gene1',
                           'init_1' = c(0,1,1),
                           'dep_1' = c(0,2,3)), row.names = F,
                file = 'testCountOne.txt')
  }
  if (!'testMasterOne.txt' %in% dir()) {
    write.table(data.table('sgrna' = paste0('g', 1:3),
                           'count' = c(0,300,300)),
                file = 'testMasterOne.txt', row.names = F)
  }
  capture.output(testDataObjAll <- suppressWarnings(DataObj$new(masterFiles = c('testMasterOne.txt'),
                                                                countFile = 'testCountOne.txt')),
                 file = 'DataObj_debugger.txt')
  # test resulting DataObj functional.
  capture.output(testModelObj <- ModelObj$new(testDataObjAll),
                 file = 'ModelObj_debugger.txt')
  checkModel <- debugModel(testModelObj, testDataObjAll, printAll = F)
  expect_true(is.data.table(testModelObj$master_freq_dt))
  expect_equal(length(testModelObj$init_scaling), 1)
  expect_false(any(is.na(c(testModelObj$init_scaling,
                           testModelObj$dep_scaling))))
  expect_true(all(sapply(c(testModelObj$init_scaling,
                           testModelObj$dep_scaling), length) == 1))

  expect_output(object = optimizeModelParameters(testDataObjAll,
                                                 testModelObj))
  testResObj <- optimizeModelParameters(testDataObjAll, testModelObj)
  expect_equal(F, is.na(testResObj$gene_results$fit_gene_param))
})

test_that('1 guide, 0 counts, dep read counts only should give DataObj warning then ModelObj error.', {
  if (!'testZeroCountsOneGuide.txt' %in% dir()) {
    write.table(data.table('guide' = paste0('g', 1),
                           'gene' = c('gene1'),
                           'dep_1' = c(0),
                           'dep_2' = c(0)),
                row.names = F,
                file = 'testZeroCountsOneGuide.txt')
  }
  if (!'testMaster.txt' %in% dir()) {
    write.table(data.table('sgrna' = paste0('g', 1:3),
                           'count' = c(200,300,300)),
                file = 'testMaster.txt', row.names = F)
  }
  expect_warning(testDataObj <- DataObj$new(masterFiles = c('testMaster.txt'),
                                            countFile = 'testZeroCountsOneGuide.txt',
                                            hasInitSeq = F),
                 'Some depleted samples have no counts in any sgRNA, check input data.')
  expect_error(testModelObj <- ModelObj$new(testDataObj),
                 'No counts in any depleted samples, provide data for negative controls.')
})
