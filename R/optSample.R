#' Function optSample
#' 
#' Given guide and gene parameters, use MLE to obtain sample parameters.
#' Sample parameters assumed to be independent, so 1D optimization used.
#' All samples present in user_DataObj are used.
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @param sample_effects Sample effects to use on initiailization.
#' @param gene_effects  By-gene essentiality parameters.
#' @param guide_efficiency By-guide efficiency parameters.
#' @param use_genes Vector of gene names to use when estimating sample parameters.
#' @param write_log Function to write to central timestamped log file.
#' @param user_DataObj DataObject with experimental data to analyze.
#' @param user_ModelObj ModelObject with preprocessed model parameters.
optSample <- function(sample_effects,
                      gene_effects,
                      guide_efficiency,
                      use_genes,
                      write_log,
                      user_DataObj,
                      user_ModelObj) {
  
  # Check for empty input.
  if (length(sample_effects)==0 | any(is.na(sample_effects))) {
    write_log('Error: Initial sample effects submitted to optSample is empty or NA')
    write_log(sample_effects)
    stop('Initial sample effects submitted to optSample is empty or NA')
  }
  if (length(use_genes)==0 | any(is.na(use_genes))) {
    write_log('Error: Gene list submitted to optSample is empty or NA')
    write_log('---')
    write_log(head(use_genes))
    write_log('---')
    stop('Error: Gene list submitted to optSample is empty or NA')
  }
  if (length(sample_effects) != ncol(user_DataObj$dep_counts)) {
    write_log('Error: length of sample effects and DataObj samples do not match.')
    stop('Length of sample effects and number of DataObj samples do not match.')
  }
  
  optObjList <- list()
  progressBar <- txtProgressBar(min = 0,
                                max = length(sample_effects),
                                initial = 0)
  for (sample_idx in seq_along(sample_effects)) {
    setTxtProgressBar(progressBar, sample_idx)
    write_log(c('Optimizing parameter for sample: ', sample_idx))
    optimOut <- capture.output(
      optObjList[[sample_idx]] <- optim(par = sample_effects[sample_idx],
                                        fn = callGetLLBySample,
                                        gene_effects = gene_effects,
                                        guide_efficiency = guide_efficiency,
                                        use_genes = use_genes,
                                        use_sample = sample_idx,
                                        user_DataObj = user_DataObj,
                                        user_ModelObj = user_ModelObj,
                                        write_log = write_log,
                                        lower = 1e-20,
                                        upper = 20,
                                        method = 'Brent',
                                        control = list(fnscale = -1,
                                                       trace = 6))
    )
    # write optim output for first parameter for debugging.
    if (sample_idx == 1) {
      write_log('=================== Sample optimization Log (iter 1) ======================')
      write_log(optimOut)
    }
  }
  return(list('optObjList' = optObjList,
              'sample_effects'= sapply(optObjList, function(i) i$par),
              'guide_efficiency' = guide_efficiency,
              'gene_effects' = gene_effects))
}