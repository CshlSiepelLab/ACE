#' Function callGetLLBySample
#' 
#' Calls Rcpp functions to evaluate the likelihood of a single SAMPLE parameter
#' given gene and guide parameters. Assumes sample independence.
#'   
#' If there is no initial sequencing, there is only lambda' (ModelObj$dep_scaling).
#' If there is initial sequencing, there is also lambda (ModelObj$init_scaling).
#' There is no scaling term for either form of the nsg prior (uniform or 
#' masterlibrary based); numerical evaluation is presumed to adjust lambda
#' to re-scale nsg to counts, assuming a linear relationship between observed counts
#' and the estimated prior as described by lambda.
#' @import data.table
#' @param sample_effect Tuple sample effect parameter.
#' @param gene_effects Numeric vector of gene essentiality parameters.
#' @param guide_efficiency Numeric vector of guide efficiency weights.
#' @param use_genes Subset of gene names to use (optional).
#' @param use_sample Sample index corresponding to sample_effect.
#' @param user_DataObj DataObject with experimental data to analyze.
#' @param user_ModelObj ModelObject with preprocessed model parameters.
# @param write_log Function to output messages to time-stamped file.
callGetLLBySample <- function(sample_effect,
                              gene_effects,
                              use_genes = NA,
                              use_sample,
                              guide_efficiency,
                              user_DataObj,
                              user_ModelObj) {
                              # write_log) {
  print('Called callGetLLBySample')
  if (any(sample_effect < 0)) sample_effect <- 0 # Prevent optim bug ignoring boundaries.
  
  # evaluate likelihood over relevant data subset.
  if (is.vector(use_genes) & !any(is.na(use_genes))) {
    useGuideCounts <- which(user_DataObj$guide2gene_map$gene %in% use_genes)
  } else if (any(is.na(use_genes))) {
    useGuideCounts <- 1:nrow(user_DataObj$guide2gene_map)
  } else {
    stop('Error: invalid gene subset argument.')
  }
  if (any(sapply(list(useGuideCounts, use_sample), length) < 1)) {
    print('Error: Samples or guides not found.')
    print(useGuideCounts[1:10])
    print('at sample index')
    print(use_sample)
    stop('Samples or guides not found.')
  }
  
  # subset count data.
  useGuides <- user_DataObj$guide2gene_map$sgrna[useGuideCounts]
  hasInit <- !any(is.na(user_ModelObj$init_scaling))
  if (hasInit) {
    gene_init_counts <- as.matrix(user_DataObj$init_counts[useGuideCounts,
                                                           (use_sample),
                                                           with=F])
    init_scaling <- sample_effect[1]
  } else {
    gene_init_counts <- NA
    init_scaling <- NA
  }
  gene_dep_counts <- as.matrix(user_DataObj$dep_counts[useGuideCounts,
                                                       (use_sample),
                                                       with=F])
  dep_sample_name <- names(user_DataObj$dep_counts)[use_sample]
  subset_dep_scaling <- user_ModelObj$dep_scaling[use_sample]
  
  # subset applicable count priors.
  masterlib <- NULL # local definition for auto-testing.
  hasGuidePrior <- user_ModelObj$guidePrior=='master_library'
  if (hasGuidePrior) {
    if (ncol(user_ModelObj$master_freq_dt)==2) { # one masterlib for all.
      masterlib_key <- 0
      gene_master_freq <- as.matrix(user_ModelObj$master_freq_dt[useGuides,2,
                                                                 with=F])
    } else {
      use_master_samples <- user_DataObj$sample_masterlib[sample %in% dep_sample_name,
                                                          unique(masterlib)]
      masterlib_key <- unlist(sapply(user_DataObj$sample_masterlib[sample %in% dep_sample_name,
                                                                   masterlib],
                                     function(i) which(use_master_samples %in% i))) - 1 # base 0 in Cpp.
      if (length(use_master_samples) < 1) {
        print('Masterlib not found')
        print(dep_sample_name)
        print(use_master_samples[1:10])
        print(user_DataObj$sample_masterlib$sample[1:10])
        stop('Could not identify paired master library.')
      }
      gene_master_freq <- as.matrix(user_ModelObj$master_freq_dt[useGuides,
                                                                 (use_master_samples),
                                                                 with=F])
    }
  } else {
    gene_master_freq <- NA
  }
  
  # subset cells infected in this sample.
  subset_cells_infected <- user_DataObj$cells_infected[use_sample]
  
  # Subset guide & gene parameter arguments; only single sample parameter being optimized.
  if (is.na(guide_efficiency)) {
    gene_guide_efficiency <- rep(1, length(useGuideCounts))
  } else {
    gene_guide_efficiency <- getGuideEfficiency(
      guide_matrix = as.matrix(user_ModelObj$guide_features[useGuideCounts,.SD]),
      feature_weight = guide_efficiency)
  }
  if (!all(use_genes %in% names(gene_effects))) {
    stop('Given gene essentiality does not contain all in use_gene.')
  } else {
    guide_essentiality <- sapply(use_genes, function(g) {
      rep(gene_effects[[g]], times = sum(user_DataObj$guide2gene_map$gene %in% g))
    })
  }
  if (hasInit & hasGuidePrior) {
    argList <- list(guide_essentiality,
                    gene_guide_efficiency,
                    sample_effect,
                    gene_init_counts,
                    gene_dep_counts,
                    user_ModelObj$mean_var_model,
                    gene_master_freq,
                    masterlib_key,
                    subset_cells_infected,
                    init_scaling, subset_dep_scaling,
                    user_ModelObj$unobserved_infected_cell_values,
                    unlist(user_ModelObj$mean_var_model_params),
                    user_ModelObj$stepSize)
    argNames <- c('guide_essentiality',
                  'gene_guide_efficiency',
                  'sample_effect',
                  'subset_init_counts',
                  'gene_dep_counts',
                  'user_ModelObj$mean_var_model',
                  'gene_master_freq',
                  'masterlib_key',
                  'subset_cells_infected',
                  'init_scaling', 'subset_dep_scaling',
                  'user_ModelObj$unobserved_infected_cell_values',
                  'unlist(user_ModelObj$mean_var_model_params)',
                  'user_ModelObj$stepSize')
    # debugging inputs:
    # debugArgList(argList, argNames)
    ll <- do.call(getLL, argList)
    
  } else if (hasInit & !hasGuidePrior) {
    noMasterArgList <- list(guide_essentiality,
                            gene_guide_efficiency,
                            sample_effect,
                            gene_init_counts,
                            gene_dep_counts,
                            user_ModelObj$mean_var_model,
                            init_scaling, subset_dep_scaling,
                            user_ModelObj$unobserved_infected_cell_values,
                            unlist(user_ModelObj$mean_var_model_params),
                            user_ModelObj$stepSize)
    argNames <- c('guide_essentiality',
                  'gene_guide_efficiency',
                  'sample_effect',
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
    noInitArgList <- list(guide_essentiality,
                          gene_guide_efficiency,
                          sample_effect,
                          gene_dep_counts,
                          user_ModelObj$mean_var_model,
                          gene_master_freq,
                          masterlib_key,
                          subset_cells_infected,
                          subset_dep_scaling,
                          user_ModelObj$unobserved_infected_cell_values,
                          unlist(user_ModelObj$mean_var_model_params),
                          user_ModelObj$stepSize)
    argNames <- c('guide_essentiality',
                  'gene_guide_efficiency',
                  'sample_effect',
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
    print('NA or inf ll returned from callGetLLBySample!')
    return(-1e20)
  }
  return(ll)
}