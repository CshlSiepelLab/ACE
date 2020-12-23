#' Function optGuide
#' 
#' Given gene and sample parameters, use MLE to obtain guide parameters.
#' @param guide_features Weights for guide efficiency features to use on initialization.
#' @param gene_effects By-gene essentiality parameters.
#' @param sample_effects By-sample parameters.
#' @param use_genes Vector of gene names to 
#' @param write_log Function to write to central timestamped log file.
#' @param user_DataObj DataObj with experimental data to analyze.
#' @param user_ModelObj ModelObj with preprocessed model parameters.
optGuide <- function(guide_features,
                     gene_effects,
                     sample_effects,
                     use_genes,
                     write_log,
                     user_DataObj,
                     user_ModelObj) {
 # Check for empty input.
 if (length(guide_features)==0 | any(is.na(guide_features))) {
   write_log('Error: Initial guide effects submitted to optGuide empty or NA')
   write_log(head(guide_features))
   stop('Initial guide effects submitted to optGuide empty or NA')
 }
  optObjList <- list()
  progressBar <- txtProgressBar(min = 0,
                                max = length(guide_features),
                                initial = 0)
  for (eg_idx in seq_along(guide_features)) {
    setTxtProgressBar(progressBar, eg_idx)
    write_log(c('Optimizing parameter for guide feature: ', eg_idx))
    optimOut <- capture.output(
      optObjList[[eg_idx]] <- optim(par = guide_features[eg_idx],
                                    fn = callGetLLByGuideFeature,
                                    user_DataObj = user_DataObj,
                                    user_ModelObj = user_ModelObj,
                                    write_log = write_log,
                                    method = 'Brent',
                                    control = list(fnscale = -1,
                                                   trace = 6))
    )
    # write optim output for first feature for debugging.
    if (sample_idx == 1) {
      write_log('=================== Sample optimization Log (iter 1) ======================')
      write_log(optimOut)
    }
  }
  return(list('optObjList'=optObjList,
              'sample_effects' = sample_effects,
              'guide_efficiency' = sapply(optObjList, function(i) i$par),
              'gene_effects' = gene_effects))
}