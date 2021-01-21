#' Function callGetLLByGene
#'
#' Calls Rcpp functions to evaluate the likelihood of a single GENE parameter
#' given sample and guide parameters.  Assumes gene independence.
#' @import data.table
#' @param geneEss Double gene essentiality parameter.
#' @param useGene String; name of gene to evaluate under given model parameters.
#' @param useSamples Numeric vector; column indices of samples to evaluate.
#' @param sample_effects Numeric vector of sample effects.
#' @param guide_efficiency Numeric vector of guide efficiency parameters.
#' @param user_DataObj DataObj class;
#' @param user_ModelObj ModelObj class.
#' @return returns log likelihood of given parameters, given data in parent ModelObjClass.
#' @export
callGetLLByGene <- function(geneEss, useGene, useSamples, sample_effects,
                            guide_efficiency,
                            user_DataObj, user_ModelObj){
  print('called callGetLLByGene')
  if (geneEss < 0) geneEss <- 0 # Boundaries set from 0 to infinity, but optim oversteps.

  # get order of rows corresponding to gene in count data.
  useGuideCounts <- which(user_DataObj$guide2gene_map$gene == useGene)
  if (any(sapply(list(useGuideCounts, useSamples), length) < 1)) {
    print('samples or guides not found.')
    print(useGuideCounts[1:10])
    print('useSamples')
    print(useSamples)
    stop('error in optimizeModelParams: genes or samples not found.')
  }
  useGuides <- user_DataObj$guide2gene_map$sgrna[useGuideCounts]
  hasInit <- is.data.table(user_DataObj$init_counts)
  hasGuidePrior <- is.data.table(user_ModelObj$master_freq_dt)
  # If hasInit design.
  if (hasInit) {
    gene_init_counts <- as.matrix(user_DataObj$init_counts[useGuideCounts,
                                                           (useSamples),
                                                           with=F])
    subset_init_scaling <- user_ModelObj$init_scaling[useSamples]
  }
  gene_dep_counts <- as.matrix(user_DataObj$dep_counts[useGuideCounts,(useSamples), with=F])
  depSampleNames <- names(user_DataObj$dep_counts)[useSamples]

  # If hasGuidePrior, model with poisson distribution of by-guide prior in master_freq.
  # Give masterlibrary frequencies subset to relevant libraries,
  # with masterlib_key of length samples indicating which idx of master_freq to use,
  # corresponding to the order of the samples in dep_counts.
  # sgrna idx column of masterlib not submitted
  
  # locally define for auto-testing.
  masterlib <- NULL 
  sgrna <- NULL
  
  if (hasGuidePrior) {
    if (ncol(user_ModelObj$master_freq_dt)==2) {
      masterlib_key <- rep(0, length(depSampleNames))
      gene_master_freq <- as.matrix(user_ModelObj$master_freq_dt[match(useGuides, sgrna),
                                                                 2,
                                                                 with=F])
    } else {
      useMasterSamples <- unique(user_DataObj$sample_masterlib[sample %in% depSampleNames,
                                                               masterlib])
      masterlib_key <- sapply(depSampleNames, function(s) {
        which(useMasterSamples == user_DataObj$sample_masterlib[sample==s,
                                                                masterlib]) - 1 # 0-based Cpp.
      })
     
      if (length(useMasterSamples) < 1) {
        print('masterlib not found')
        print(useMasterSamples[1:10])
        print(user_DataObj$sample_masterlib[1:10])
        stop('Could not identify paired master library.')
      }
      if (useMasterSamples[masterlib_key[1]+1] != 
          user_DataObj$sample_masterlib[sample==depSampleNames[1], masterlib]){
        stop('Incorrectly subset master library.')
      }

      gene_master_freq <- as.matrix(user_ModelObj$master_freq_dt[match(useGuides, sgrna),
                                                                 (useMasterSamples),
                                                                 with=F])
    }
    subset_cells_infected <- user_DataObj$cells_infected[useSamples]
  }

  # Convert guide efficiency feature weights to % effective by guide.
  if (is.na(guide_efficiency)) {
    gene_guide_efficiency <- rep(1, length(useGuideCounts))
  } else if (is.vector(guide_efficiency)) {
    gene_guide_efficiency <- getGuideEfficiency(
      guide_matrix = as.matrix(user_ModelObj$guide_features[useGuideCounts,.SD]),
      feature_weight = guide_efficiency)
  } 
  
  # Format arguments for cpp code - note 0-based indexing!!
  gene_essentiality <- rep(geneEss, length(gene_guide_efficiency))
  subset_sample_effects <- sample_effects[useSamples]
  subset_dep_scaling <- user_ModelObj$dep_scaling[useSamples]
  if (hasInit & hasGuidePrior) {
    argList <- list(gene_essentiality,
                    gene_guide_efficiency,
                    subset_sample_effects,
                    gene_init_counts,
                    gene_dep_counts,
                    user_ModelObj$mean_var_model,
                    gene_master_freq,
                    masterlib_key,
                    subset_cells_infected,
                    subset_init_scaling, subset_dep_scaling,
                    user_ModelObj$unobserved_infected_cell_values,
                    unlist(user_ModelObj$mean_var_model_params),
                    user_ModelObj$stepSize)
    argNames <- c('gene_essentiality',
                  'gene_guide_efficiency',
                  'subset_sample_effects',
                  'gene_init_counts',
                  'gene_dep_counts',
                  'user_ModelObj$mean_var_model',
                  'gene_master_freq',
                  'masterlib_key',
                  'subset_cells_infected',
                  'subset_init_scaling', 'subset_dep_scaling',
                  'user_ModelObj$unobserved_infected_cell_values',
                  'unlist(user_ModelObj$mean_var_model_params)',
                  'user_ModelObj$stepSize')
    # debugging inputs:
    checkArgList <- argList
    # debugArgList(argList, argNames)
    ll <- do.call(getLL, argList)

  } else if (hasInit & !hasGuidePrior) {
    noMasterArgList <- list(gene_essentiality,
                            gene_guide_efficiency,
                            subset_sample_effects,
                            gene_init_counts,
                            gene_dep_counts,
                            user_ModelObj$mean_var_model,
                            subset_init_scaling, subset_dep_scaling,
                            user_ModelObj$unobserved_infected_cell_values,
                            unlist(user_ModelObj$mean_var_model_params),
                            user_ModelObj$stepSize)
    argNames <- c('gene_essentiality',
                  'gene_guide_efficiency',
                  'subset_sample_effects',
                  'gene_init_counts',
                  'gene_dep_counts',
                  'user_ModelObj$mean_var_model',
                  'subset_init_scaling', 'subset_dep_scaling',
                  'user_ModelObj$unobserved_infected_cell_values',
                  'unlist(user_ModelObj$mean_var_model_params)',
                  'user_ModelObj$stepSize')
     # debugArgList(noMasterArgList, argNames)
    checkArgList <- noMasterArgList
    ll <- do.call(getLLNoMaster, noMasterArgList)

  } else if (hasGuidePrior & !hasInit) {
    noInitArgList <- list(gene_essentiality,
                          gene_guide_efficiency,
                          subset_sample_effects,
                          gene_dep_counts,
                          user_ModelObj$mean_var_model,
                          gene_master_freq,
                          masterlib_key,
                          subset_cells_infected,
                          subset_dep_scaling,
                          user_ModelObj$unobserved_infected_cell_values,
                          unlist(user_ModelObj$mean_var_model_params),
                          user_ModelObj$stepSize)
    argNames <- c('gene_essentiality',
                  'gene_guide_efficiency',
                  'subset_sample_effects',
                  'gene_dep_counts',
                  'user_ModelObj$mean_var_model',
                  'gene_master_freq',
                  'masterlib_key',
                  'subset_cells_infected',
                  'subset_dep_scaling',
                  'user_ModelObj$unobserved_infected_cell_values',
                  'unlist(user_ModelObj$mean_var_model_params)',
                  'user_ModelObj$stepSize')
    # debugArgList(noInitArgList, argNames)
    checkArgList <- noInitArgList
    ll <- do.call(getLLNoInit, noInitArgList)
  } else {
    stop('invalid experimental design. Must have master/init and depleted sets.')
  }
  if (is.na(ll) | is.infinite(ll)) {
    debugArgList(checkArgList, argNames, printArgList = T)
    print('NA or inf ll returned from callGetLLByGene!')
    return(-1e20)
  }
  return(ll)
}
