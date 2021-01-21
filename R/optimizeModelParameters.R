#' Function optimizeModelParameters
#'
#' Given model object and data object, create a results object.
#' Mainly runs optim over sample (sub)set(s).
#' in final version, have as part of R6 object which inherits from ModelObj and DataObj;
#' for debugging, run as stand-alone with R6 as input.
#' @import data.table
#' @import stringr
#' @importFrom stats optim
#' @importFrom utils capture.output head timestamp write.table
#' @importFrom Rcpp sourceCpp
## usethis namespace: start
#' @useDynLib ACER, .registration = TRUE
## usethis namespace: end
#' @param user_DataObj DataObj with all experimental data.
#' @param user_ModelObj ModelObj with simple inferred parameters for model.
#' @param fit_sample_param User option to fit by-sample parameters.
#' @param fit_guide_param User option to fit by-guide feature weights.
#' @param subset_genes optional list of start & end points of gene indices to subset.
#' @param converge_ll Tolerated variation in likelihood to end optimization; default 100.
#' @param max_iter Maximum number of optimization iterations, regardless of convergence. Default 10.
#' @param ncpus Number of threads to use in parallel optimization; auto-detected if not specified.
#' @return Data.table with Gene name, optimized essentiality,  guide_efficiency, sample_response,
#'          guide_covariates, gene_essentiality likelihood, gene_ess CI (hessian),; one DT each for
#'          "all", "test", and "control" conditions, and an optional 4th for differential depletion
#'          effect sizes.
#' @export

optimizeModelParameters <- function(user_DataObj, user_ModelObj, 
                                    fit_sample_param = F,
                                    fit_guide_param = F,
                                    subset_genes = NA,
                                    converge_ll = 100,
                                    max_iter = 10,
                                    ncpus = NA) {
  # Set local definitions to prevent R check note due to data.table syntax.
  test_fit <- NULL
  ctrl_fit <- NULL
  all_fit <- NULL
  gene <- NULL
  fit_gene_param <- NULL
  
  # -------------------------------- Functions ------------------------------------- #
  
  write_log_set_file <- function(message_vector, log_file) {
    if (is.atomic(message_vector)) {
      cat(message_vector,'\n', file=log_file, append=T)
    } else {
      suppressWarnings(write.table(message_vector, file = log_file, 
                                   col.names=T, append=T))
    }
  }

  # -------------------------------- Main Method -------------------------------------

  # Create log file. 
  if (!dir.exists('ACE_output_data')) dir.create('ACE_output_data')
  tStamp <- paste(unlist(str_split(Sys.time(), ' |:')), collapse='_')
  log_file <- file(file.path('ACE_output_data',paste0('ACE_optim_log_', tStamp,
                                                      '.txt')),
                   open='w+')
  on.exit(close(log_file), add=T)
  write_log <- function(message_vector) {
    write_log_set_file(message_vector, log_file)
    return()
  }
  
  # Check optimization options valid.
  if (fit_guide_param) {
    if (!is.data.table(user_ModelObj$guide_features)) {
      fit_guide_param <- F
      message('Unable to fit guide parameters, no features contained in ModelObj.')
      write_log('Unable to fit guide parameters, no features contained in ModelObj.')
    } else {
      message('Fitting guide parameters.')
      write_log('Fitting guide parameters.')
    }
  }
  
  # Determine indices of genes to optimize. 
  geneList <- unique(user_DataObj$guide2gene_map$gene)
  numGenes <- length(geneList)
  if (is.list(subset_genes) & length(subset_genes)==2) {
    if (any(c(subset_genes[[1]] < 1,
              subset_genes[[2]] > nrow(user_DataObj$dep_counts),
              length(subset_genes) > 2))) {
      message('Only ',numGenes, ' genes queried.')
      write_log(c('Only ',numGenes, ' genes queried.'))
      write_log('Error: invalid gene subset specified.')
      stop('invalid gene subset specified.')
    }
    subsetG <- do.call(seq, subset_genes)
  } else if (is.vector(subset_genes) & all(!is.na(subset_genes))) {
    if (is.numeric(subset_genes)) {
      subsetG <- subset_genes
    } else {
      subsetG <- which(geneList %in% subset_genes)
    }
  } else if (is.na(subset_genes[1])) {
    subsetG <- 1:numGenes
  } else {
    write_log('Error: Invalid gene subset argument:')
    write_log(subset_genes)
    stop('invalid gene subset argument.')
  }
  # Initialize guide & gene parameters for first round of gene essentiality 
  # estimation as perfectly efficient.
  guide_efficiency <- NA
  guide_essentiality <- rep(1, nrow(user_DataObj$dep_counts))
  sampleSubsets <- list("all"=1:ncol(user_DataObj$dep_counts))
  startEss <- rep(1, length(geneList))
  names(startEss) <- geneList
  gene_essentiality <- startEss
  sample_effects <- rep(1, ncol(user_DataObj$dep_counts))
  if (is.vector(user_ModelObj$test_sample_subtype_cols, mode='integer')) {
    sampleSubsets$test <- user_ModelObj$test_sample_subtype_cols
    sampleSubsets$ctrl <- sampleSubsets$all[!sampleSubsets$all %in% sampleSubsets$test]
  }

  # 3 - optimize by-gene params;
  # Determine differential depletion in sample subsets by fitting TWO essential
  # parameters per gene, one global and one subset-specific.
  optFgObj <- optFg(startEss, sampleSubsets, sample_effects,
                    guide_efficiency, geneList[subsetG], write_log,
                    user_DataObj = user_DataObj,
                    user_ModelObj = user_ModelObj,
                    ncpus = ncpus)
  message('Initial by-gene parameter optimization complete')
  optObjList <- optFgObj[[1]]
  nullLogLikeList <- optFgObj[[2]] 
  
  # 4 - iteratively optimize sample-, guide-, and gene-specific parameters.
  # Optional.
  iterNum <- 0
  if (fit_sample_param | fit_guide_param) {
    continueOptimizing <- T
    isConverged <- F
    while (continueOptimizing) {
      old_likelihood_estimate <- sum(optObjList$gene_results$ll)
      
      if (fit_sample_param) {
        optSampleObj <- optSample(optFgObj$sample_effects,
                                  optFgObj$gene_effects,
                                  optFgObj$guide_efficiency,
                                  geneList[subsetG],
                                  write_log,
                                  user_DataObj = user_DataObj,
                                  user_ModelObj = user_ModelObj)
      } else {
        optSampleObj <- optFgObj # just a pointer, not a clone.
      }
      
      if (fit_guide_param) {
        optGuideObj <- optGuide(optFgObj$guide_efficiency,
                                optFgObj$gene_effects,
                                optFgObj$sample_effects,
                                geneList[subsetG],
                                write_log,
                                user_DataObj = user_DataObj,
                                user_ModelObj = user_ModelObj)
      } else {
        optGuideObj <- optSampleObj # just a pointer, not a clone.
      }
      # update gene efficiency estimates based on sample and guide effects,
      # for all sample subsets.
      optFgObj <- optFg(optFgObj$gene_effects, 
                        sampleSubsets, 
                        optSampleObj$sample_effects,
                        optGuideObj$guide_efficiency, 
                        geneList[subsetG], write_log,
                        user_DataObj = user_DataObj,
                        user_ModelObj = user_ModelObj,
                        ncpus = ncpus)
      iterNum <- iterNum + 1
      
      # Exit conditions.
      delta_ll <- abs(old_likelihood_estimate - sum(optObjList$gene_results$ll))
      
      if (iterNum == max_iter | delta_ll < converge_ll) {
        write_log('Exit conditions met, converged: ', delta_ll < converge_ll)
        continueOptimizing <- F
      }
    }
  }      
  message(paste0('Optimization complete, iteration: ', iterNum))
  write_log(paste0('Optimization complete, iteration: ', iterNum))
  optObjList <- optFgObj[[1]]
  nullLogLikeList <- optFgObj[[2]]
  gene_essentiality[subsetG] <- sapply(optObjList[['all']], function(i) i$par)
  
  # TODO: final_ReturnObj <- ReturnObj$new(optObjList, user_DataObj, user_ModelObj)
  # @return List with:
  #         -- gene_results
  #         -- diff_genes
  #         Data.table with Gene name, optimized essentiality,  guide_efficiency, sample_response,
  #          guide_covariates, gene_essentiality likelihood, gene_ess CI (hessian),; one DT each for
  #          "all", "test", and "control" conditions, and an optional 4th for differential depletion
  #          effect sizes.
  # return 'gene results' as an all-sample fit.
  geneList <- geneList[geneList %in% names(optObjList[['all']])]
  gene_results <- data.table("gene" = geneList,
                             "fit_gene_param" = sapply(geneList, function(i)
                               optObjList[['all']][[i]]$par),
                             "log_ll_fit" = sapply(geneList, function(i)
                               optObjList[['all']][[i]]$value),
                             "log_ll_null" = sapply(geneList, function(i)
                               nullLogLikeList[['all']][[i]]))
  fit_hess <- lapply(geneList, function(g) optObjList[['all']][[g]]$hessian)
  try(gene_results[, 'fit_se_all' := sapply(seq_along(fit_hess),
                                            function(i) get95Int(fit_hess[[i]],
                                                                 gene_results$fit_gene_param[i])[[3]])])
  write_log('Sample results:')
  write_log(head(gene_results))

  ## Get results for differential essentiality between test & ctrl samples. ##
  diff_genes <- NA
  # log_ll_null = log_ll_all; log_ll_alt = log_ll_test + log_ll_ctrl; df=1
  if (length(sampleSubsets) > 1) {
    diff_genes <- data.table("gene" = geneList,
                             "test_fit" = sapply(geneList, function(i)
                                                 optObjList[['test']][[i]]$par),
                             "test_ll" = sapply(geneList, function(i)
                                                optObjList[['test']][[i]]$value),
                             "ctrl_fit" = sapply(geneList, function(i)
                                                 optObjList[['ctrl']][[i]]$par),
                             "ctrl_ll" = sapply(geneList, function(i)
                                                optObjList[['ctrl']][[i]]$value),
                             "all_fit" = gene_results$fit_gene_param,
                             "all_ll" = gene_results$log_ll_fit,
                             "test_null_ll" = sapply(geneList, function(i)
                               nullLogLikeList[['test']][[i]]),
                             "ctrl_null_ll" = sapply(geneList, function(i)
                               nullLogLikeList[['ctrl']][[i]]))
    diff_genes[, "diff_depletion" := test_fit/ctrl_fit]
    dep_hess <- lapply(geneList, function(g) optObjList[['test']][[g]]$hessian)
    ctrl_hess <- lapply(geneList, function(g) optObjList[['ctrl']][[g]]$hessian)
    
    # scale absolute essentiality value according to negative controls, if given.
    if (user_ModelObj$use_neg_ctrl & any(user_ModelObj$neg_ctrls %in% 
                                         diff_genes$gene)) {
      all_scale <- diff_genes[gene %in% user_ModelObj$neg_ctrls, median(all_fit)]
      test_scale <- diff_genes[gene %in% user_ModelObj$neg_ctrls, median(test_fit)]
      ctrl_scale <- diff_genes[gene %in% user_ModelObj$neg_ctrls, median(ctrl_fit)]
      gene_results[, fit_gene_param := fit_gene_param/all_scale]
      write_log('Inferred sample parameters:')
      write_log(c(all_scale, test_scale, ctrl_scale))
      diff_genes[, all_fit := all_fit/all_scale]
      diff_genes[, test_fit := test_fit/test_scale]
      diff_genes[, ctrl_fit := ctrl_fit/ctrl_scale]
    }
  } else{
    dep_hess <- NA
    ctrl_hess <- NA
    diff_genes <- NA
  }

  # TODO: return by-guide everything.
  guide_results <- data.table("guide_feature",
                              "fit_guide_param",
                              "fit_CI_upper",
                              "fit_CI_lower")
  # Return by-sample weights used.
  sample_results <- data.table("sample" = sampleSubsets[['all']],
                               "fit_sample_response" = sample_effects)

  readMe <- paste0("log_ll are the log of the maximum values produced by callGetLLByGene,",
                "with null from no gene effect.")
  final_ReturnObj <- list("gene_results"=gene_results,
                          "guide_results"=guide_results,
                          "sample_results"=sample_results,
                          "diff_genes"=diff_genes,
                          "readMe"=readMe,
                          'fit_hess' = fit_hess,
                          'dep_hess'=dep_hess,
                          'ctrl_hess' = ctrl_hess)
  return(final_ReturnObj)
}

# -------------------------------- End of Script -------------------------------------
