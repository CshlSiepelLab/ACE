context("Parameter Optimization")
library(ACER)
setup({
  write.table(data.table('guide' = paste0('g', 1),
                         'gene' = c('gene1'),
                         'dep_1' = c(0),
                         'dep_2' = c(0)),
              row.names = F,
              file = 'testZeroCountsOneGuide.txt')
  write.table(data.table('guide' = paste0('g', 1:3),
                         'gene' = 'gene1',
                         'dep_1' = c(0,2,3),
                         'dep_2' = c(0,4,5)),
              row.names = F,
              file = 'testZeroCountsInOne.txt')
  write.table(data.table('sgrna' = paste0('g', 1:3),
                         'count' = c(200,300,300)),
              file = 'testMaster.txt', row.names = F)
  write.table(data.table('sgrna' = paste0('g', 1:3),
                         'count' = c(2000,3000,3000)),
              file = 'testMaster2.txt', row.names = F)
})
teardown({
  file.remove('testZeroCountsOneGuide.txt', 'testZeroCountsInOne.txt',
              'testMaster.txt', 'testMaster2.txt')
  file.remove(sapply(dir('ACE_output_data'),function(f) {
    file.path('ACE_output_data', f)}), 'ACE_output_data')
})

# For testing external function interfaces.
test_that("get_pois_muVec returns same results as R::dpois", {
  expect_equal(get_pois_muVec(0, 1:1e4), dpois(0, 1:1e4, log=T))
  expect_equal(get_pois_muVec(0, seq(1, 2, .01)), dpois(0, seq(1, 2, .01), log=T))
  expect_equal(get_pois_muVec(0, 0), dpois(0,0, log=T))
})

test_that('Call getLLNoInit with zero counts in all depleted samples', {
  argList <- list('gene_essentiality' = 0,
                  'gene_guide_efficiency' = 1,
                  'subset_sample_effects' = rep(1, 10),
                  'gene_dep_counts' = matrix(rep(0, 10), nrow=1),
                  'user_ModelObj$mean_var_model' = 1,
                  'gene_master_freq' = matrix(c(-17.0605, -17.22331),
                                                 nrow=1, ncol=2),
                  'masterlib_key' = c(0, 0, 1,1,1, 0,0,0,1,1),
                  'subset_cells_infected' = rep(90710000, 10),
                  'subset_dep_scaling' = rep(10, 10),
                  'user_ModelObj$unobserved_infected_cell_values' = seq(1, 1e4, 10),
                  'unlist(user_ModelObj$mean_var_model_params)' = c(0,0),
                  'stepSize' = 10)
  names(argList) <- NULL
  expect_equal(2.19, round(do.call(getLLNoInit, argList), digits = 2))
  expect_false(is.na(do.call(getLLNoInit, argList)))
})

test_that('0 counts in depleted samples, no init, passed to getLLNoInit', {

  testDataObj <- DataObj$new(masterFiles = c('testMaster.txt'),
                             countFile = 'testZeroCountsInOne.txt',
                             hasInitSeq = F)
  expect_equal(3, nrow(testDataObj$dep_counts))
  testModelObj <- ModelObj$new(testDataObj)
  expect_output(testResObj <- optimizeModelParameters(testDataObj,
                                        testModelObj))
  expect_false(is.na(testResObj$gene_results$fit_gene_param))

})

test_that('Multiple masterlibraries can be included without init seq.', {


  expect_output(testDataObj <- DataObj$new(masterFiles = c('testMaster.txt',
                                                           'testMaster2.txt'),
                                           countFile = 'testZeroCountsInOne.txt',
                                           hasInitSeq = F))
})

# what happens when you try to optimize with 0 counts?
test_that('Require counts in some guides in depleted samples if no init seq.', {
  expect_warning(testDataObj <- DataObj$new(masterFiles = c('testMaster.txt'),
                                            countFile = 'testZeroCountsOneGuide.txt',
                                            hasInitSeq = F),
                 'Some depleted samples have no counts in any sgRNA, check input data.')
  expect_error(testModelObj <- ModelObj$new(testDataObj),
               'No counts in any depleted samples, provide data for negative controls.')
})

# optim() calls values outside of the parameter bounds during gradient estimation.
test_that('Negative and positive ess values should not be accepted by callGetLLByGene.', {

  # no warnings if no counts just in depleted sample.
  expect_output(testDataObj <- DataObj$new(masterFiles = c('testMaster.txt'),
                                            countFile = 'testZeroCountsInOne.txt',
                                            hasInitSeq = F))
  testModelObj <- ModelObj$new(testDataObj)
  baseParams <- list(useGene = 'gene1',
                     useSamples = 1,
                     sample_effects = 1,
                     guide_efficiency = c(1,1,1),
                     user_DataObj = testDataObj,
                     user_ModelObj = testModelObj)
  expect_output(res1 <- do.call(callGetLLByGene, append(list('geneEss' = 1), baseParams)))
  expect_output(res2 <- do.call(callGetLLByGene, append(list('geneEss' = 0), baseParams)))
  expect_output(res3 <- do.call(callGetLLByGene, append(list('geneEss' = -1), baseParams)))
  expect_output(res4 <- do.call(callGetLLByGene, append(list('geneEss' = 2), baseParams)))
  # expect values returned to be the same at and below/above the lower/upper
  # boundaries, respectively.
  expect_equal(res2, res3)
  # based on current parameterization of fit_gene_param, may change if max is 1:
  # expect_equal(res1, res4)
})
