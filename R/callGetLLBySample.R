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
#' @param write_log Function to output messages to time-stamped file.
callGetLLBySample <- function(sample_effect,
                              gene_effects,
                              use_genes = NA,
                              use_sample,
                              guide_efficiency,
                              user_DataObj,
                              user_ModelObj,
                              write_log) {
  write_log('Called callGetLLBySample')
  if (sample_effect < 0) sample_effect <- 0 # Prevent optim bug ignoring boundaries.
  
  # evaluate likelihood over relevant data subset.
  if (is.vector(use_genes)) {
    useGuideCounts <- which(user_DataObj$guide2gene_map$gene %in% useGene)
  } else if (is.na(use_genes)) {
    useGuideCounts <- 1:nrow(user_DataObj$guide2gene_map)
  } else {
    stop('Error: invalid gene subset argument.')
  }
  if (any(sapply(list(useGuideCounts, use_samples), length) < 1)) {
    write_log('Error: Samples or guides not found.')
    write_log(useGuideCounts[1:10])
    write_log('at sample index')
    write_log(use_sample)
    stop('Samples or guides not found.')
  }
  
  # subset count data.
  useGuides <- user_DataObj$guide2gene_map$sgrna[useGuideCounts]
  hasInit <- is.vector(user_ModelObj$init_scaling)
  if (hasInit) {
    gene_init_counts <- as.matrix(user_DataObj$init_counts[useGuideCounts,
                                                           (use_sample),
                                                           with=F])
    init_scaling <- sample_effect[1]
  }
  gene_dep_counts <- as.matrix(user_DataObj$dep_counts[useGuideCounts,
                                                       (use_sample),
                                                       with=F])
  dep_sample_name <- names(user_DataObj$dep_counts)[use_sample]
  
  # subset applicable count priors.
  masterlib <- NULL # local definition for auto-testing.
  if (user_ModelObj$guidePrior=='master_library') {
    if (ncol(userModelObj$master_freq_dt)==2) { # one masterlib for all.
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
        write_log('Masterlib not found')
        write_log(dep_sample_name)
        write_log(use_master_samples[1:10])
        write_log(user_DataObj$sample_masterlib$sample[1:10])
        stop('Could not identify paired master library.')
      }
      # data.table strips '.txt' from column names in ModelObj.
      use_master_samples <- sapply(use_master_samples, function(i) {
        strsplit(i, '.txt')[[1]]
        })
      gene_master_freq <- as.matrix(user_ModelObj$master_freq_dt[useGuides,
                                                                 (use_master_samples),
                                                                 with=F])
    }
    
    # subset cells infected in this sample.
    subset_cells_infected <- user_DataObj$cells_infected[use_sample]
    
  }
  
}