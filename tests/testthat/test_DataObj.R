context('Testing DataObj object.')
library(ACER)
setup({
  write.table(data.table('guide' = paste0('g', 1:3),
                         'gene' = 'gene1',
                         'init_1' = c(0,1,1),
                         'dep_1' = c(0,2,3)),
              quote = F, row.names = F,
              file = 'test_CountOne.txt')
  write.table(data.table('sgrna' = paste0('g', 1:3),
                         'count' = c(0,300,300)),
              file = 'test_MasterOne.txt',
              quote = F, row.names = F)
  write.table(data.table('sgrna' = paste0('g', 1:3),
                         'gene' = 'gene1',
                         'dep_1' = c(0,3,3),
                         'dep_2' = c(0, 50, 50)),
              row.names = F, quote = F,
              file = 'test_CountOneDep.txt')
  write.table(data.table('guide' = paste0('g', 1:3),
                         'gene' = 'gene1',
                         'init_1' = rep(0,3),
                         'dep_1' = rep(0,3)), row.names = F,
              file = 'test_CountNull.txt')
  write.table(data.table('sgrna' = paste0('g', 1:3),
                         'count' = rep(0,3)),
              file = 'test_MasterNull.txt', row.names = F)
  write.table(data.table('guide' = paste0('g', 1:3),
                         'gene' = 'gene1',
                         'init_1' = c(2,2,2),
                         'dep_1' = c(1,1,1),
                         'init2' = c(3,0,3),
                         'dep_2' = c(2,0,2)),
              file = 'test_DropSample.txt', row.names = F)
})

teardown({
  file.remove('test_CountOne.txt',
              'test_MasterOne.txt',
              'test_CountOneDep.txt',
              'test_CountNull.txt',
              'test_MasterNull.txt',
              'test_DropSample.txt')

  file.remove(sapply(dir('ACE_output_data'),function(f) {
    file.path('ACE_output_data', f)}), 'ACE_output_data')
  file.remove('DataObj_debugger.txt')
})

test_that('Test no counts in any gene in any sample throws an error.', {

  expect_error(capture_warnings(DataObj$new(masterFiles = c('test_MasterNull.txt'),
                                countFile = 'test_CountNull.txt')),
               'No valid data submitted.')
})
test_that('Blank guides should be removed from count data with a warning.', {

  warnAll <- capture_warnings(testDataObjAll <- DataObj$new(masterFiles = c('test_MasterOne.txt'),
                                              countFile = 'test_CountOne.txt'))

  warnMaster <- capture_warnings(testDataObjMaster <- DataObj$new(masterFiles = c('test_MasterOne.txt'),
                                                                  countFile = 'test_CountOneDep.txt',
                                                                  hasInitSeq = F))
  warnInit <- capture_warnings(testDataObjInit <- DataObj$new(countFile = 'test_CountOne.txt'))

  expect_equal(2, nrow(testDataObjAll$dep_counts))
  expect_equal(2, nrow(testDataObjMaster$dep_counts))
  expect_equal(2, nrow(testDataObjInit$dep_counts))

  checkData <- debugData(testDataObjAll, printAll = F)
  expect_equal(0, length(checkData$flag))
  expect_equal(length(testDataObjAll$init_total_reads), 1)
})

test_that('Test guides missing from only one sample are included.', {

  expect_output(DataObj$new(countFile = 'test_DropSample.txt'))
})
