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
#' @param write_log Function with time-stamped file to output messages.
#' @return returns log likelihood of given parameters, given data in parent ModelObjClass.
#' @export
callGetLLByGene <- function(geneEss, useGene, useSamples, sample_effects,
                            guide_efficiency,
                            user_DataObj, user_ModelObj, write_log) {
  write_log('called callGetLLByGene')
  if (geneEss < 0) geneEss <- 0

  # get order of rows corresponding to gene in count data.
  useGuideCounts <- which(user_DataObj$guide2gene_map$gene == useGene)
  if (any(sapply(list(useGuideCounts, useSamples), length) < 1)) {
    write_log('samples or guides not found.')
    write_log(useGuideCounts[1:10])
    write_log('useSamples')
    write_log(useSamples)
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
  # DON'T expand master_freq into a matrix guides x samples;
  # give masterlibrary of just useGuides and unique(masterlib[depSampleNames]),
  # with masterlib_key of length samples indicating which idx of master_freq to use.
  # sgrna idx column of masterlib not submitted
  masterlib <- NULL # locally define for auto-testing.
  if (hasGuidePrior) {
    if (ncol(user_ModelObj$master_freq_dt)==2) {
      masterlib_key <- rep(0, length(depSampleNames))
      gene_master_freq <- as.matrix(user_ModelObj$master_freq_dt[useGuides, 2,
                                                                 with=F])
    } else {
      useMasterSamples <- user_DataObj$sample_masterlib[sample %in% depSampleNames,
                                                        unique(masterlib)]
      masterlib_key <- unlist(sapply(user_DataObj$sample_masterlib[sample %in% depSampleNames,
                                                                   masterlib],
                                     function(i) which(useMasterSamples %in% i))) - 1 # base 0 in Cpp.
    if (length(useMasterSamples) < 1) {
      write_log('masterlib not found')
      write_log(useMasterSamples[1:10])
      write_log(user_DataObj$sample_masterlib[1:10])
      stop('Could not identify paired master library.')
    }
    # data.table strips '.txt' from column names somewhere in MOdelobj.
    useMasterSamples <- sapply(useMasterSamples, function(i) strsplit(i, '.txt')[[1]])
    gene_master_freq <- as.matrix(user_ModelObj$master_freq_dt[useGuides,
                                                               (useMasterSamples),
                                                               with=F])
    }
    subset_cells_infected <- user_DataObj$cells_infected[useSamples]
  }

  # Format arguments for cpp code - note 0-based indexing!!
  gene_guide_efficiency <- guide_efficiency[useGuideCounts]
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
    ll <- do.call(getLLNoInit, noInitArgList)
  } else {
    stop('invalid experimental design. Must have master/init and depleted sets.')
  }
  if (is.na(ll) | is.infinite(ll)) {
    write_log('NA or inf ll returned!')
    return(-1e20)
  }
  return(ll)
}
