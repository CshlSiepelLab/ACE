context('Test Model Object')
library(ACER)

setup({
  write.table(data.table('guide' = paste0('g', 1:3),
                         'gene' = 'gene1',
                         'init_1' = c(0,1,1),
                         'dep_1' = c(0,2,3)), row.names = F,
              file = 'testCountOne.txt')
  write.table(data.table('sgrna' = paste0('g', 1:3),
                         'count' = c(0,300,300)),
              file = 'testMasterOne.txt', row.names = F)
  write.table(data.table('guide' = paste0('g', 1),
                         'gene' = c('gene1'),
                         'dep_1' = c(0),
                         'dep_2' = c(0)),
              row.names = F,
              file = 'testZeroCountsOneGuide.txt')
  write.table(data.table('sgrna' = paste0('g', 1:3),
                         'count' = c(200,300,300)),
              file = 'testMaster.txt', row.names = F)
  write(x='Debugging Messages from processing master library input files.',
        file = 'getMasterLib_debugger.txt')
  write(x='Debugging Messages from ModelObj Class',
        file = 'ModelObj_debugger.txt')
})

teardown({
  file.remove('testCountOne.txt', 'testMasterOne.txt',
              'testZeroCountsOneGuide.txt',
              'testMaster.txt')
  file.remove(sapply(dir('ACE_output_data'),function(f) {
    file.path('ACE_output_data', f)}), 'ACE_output_data')
  file.remove('ModelObj_debugger.txt')
  file.remove('getMasterLib_debugger.txt')
})

test_that('Master library replicate handling in getMasterLib.', {
  testMasterCountList <- list('masterA' = data.table('sgrna'=paste0('g', 1:1e4),
                                                     'rep1' = rep(0:9, 1e3)) ,
                              'masterB' = data.table('sgrna' = paste0('g', 1:1e4),
                                                     'repli1' = 1:1e4,
                                                     'repli2' = (1+2e4):3e4),
                              'masterC' = data.table('sgrna' = paste0('g', 1:1e3),
                                                     'rep1' = rep(0, 1e3)))
  log_file <- file('getMasterLib_debugger.txt', open = 'w+')
  on.exit(close(log_file))
  testLog <- function(message_vector) {
    if (is.atomic(message_vector)) {
      cat(message_vector,'\n', file=log_file, append=T)
    } else {
      suppressWarnings(write.table(message_vector, file = log_file, 
                                   col.names=T, append=T))
    }
  }
  testLog(head(testMasterCountList[[1]]))
  testMasterFreq <- getMasterLib(countList = testMasterCountList, 
                                                testLog)
  noInfValuesInMasterFreq <- !testMasterFreq[, any(sapply(.SD, is.infinite)),
                                             .SDcols = -'sgrna']
  expect_true(noInfValuesInMasterFreq)
  noNaValuesInMasterFreq <- !testMasterFreq[, any(sapply(.SD, is.na)),
                                            .SDcols = -'sgrna']
  expect_true(noNaValuesInMasterFreq)
  expect_true(is.data.table(testMasterFreq))
})

test_that('ModelObj with blank guides has problem.', {

  testDataObjAll <- suppressWarnings(DataObj$new(masterFiles = c('testMasterOne.txt'),
                                                                countFile = 'testCountOne.txt'))
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

  expect_warning(testDataObj <- DataObj$new(masterFiles = c('testMaster.txt'),
                                            countFile = 'testZeroCountsOneGuide.txt',
                                            hasInitSeq = F),
                 'Some depleted samples have no counts in any sgRNA, check input data.')
  expect_error(testModelObj <- ModelObj$new(testDataObj),
                 'No counts in any depleted samples, provide data for negative controls.')
})
