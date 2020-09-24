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
#' @param fit_sample_parameter User option to fit by-sample parameters.
#' @param subset_genes optional list of start & end points of gene indices to subset.
#' @return Data.table with Gene name, optimized essentiality,  guide_efficiency, sample_response,
#'          guide_covariates, gene_essentiality likelihood, gene_ess CI (hessian),; one DT each for
#'          "all", "test", and "control" conditions, and an optional 4th for differential depletion
#'          effect sizes.
#' @export

optimizeModelParameters <- function(user_DataObj, user_ModelObj, fit_sample_parameter=F,
                                    subset_genes = NA) {
  # Set local definitions to prevent R check note due to data.table syntax.
  sample_effects <- NULL
  test_fit <- NULL
  ctrl_fit <- NULL
  all_fit <- NULL
  gene <- NULL
  fit_gene_param <- NULL
  
  # -------------------------------- Functions ------------------------------------- #
  # given guide feature matrix, apply by-feature weight (features apply across guides),
  # and transform the resulting covariate score to binary active/inactive or
  # to percent efficiency.
  getGuideEfficiency <- function(covar_weights) {
    covar_score <- as.vector(user_ModelObj$guide_covar %*% t(covar_weights))
    guide_efficiency <- ifelse(covar_score > 0.5, 1, 0)
    return(guide_efficiency)
  }

  # assume each gene is independent, and optimize over a single gene at a time.
  # while this single parameter optimization is very fast, this for loop *can* be parallelized.
  # produces list of optim objects.
  # Assume phi's positive selection impact capped at 20-fold.
  #mu_1 = lambda1[s]*n
  #mu_2 = lambda2[s]* mu_1 * (1 - eps + eps * phi)
  # (see src/getLL.cpp for the actual math).
  optFg <- function(startEss, sampleSubsets, sample_effects,
                    guide_efficiency, geneList) {
    optObjList <- list()
    nullLogLikeList <- list()
    tStamp <- paste(unlist(stringr::str_split(Sys.time(), ' ')), collapse='_')
    for (subset in seq_along(sampleSubsets)) {
      for (gene_idx in seq_along(geneList)) {
        gene <- geneList[[gene_idx]]
        cat("optimizing parameters for gene: ", gene, '\n')
        optimOut <- capture.output(
          optObjList[[names(sampleSubsets)[subset]]][[gene]] <- optim(par = startEss[gene],
                                                                      fn = callGetLLByGene,
                                                                      useGene = gene,
                                                                      useSamples = sampleSubsets[[subset]],
                                                                      sample_effects = sample_effects,
                                                                      guide_efficiency = guide_efficiency,
                                                                      user_DataObj = user_DataObj,
                                                                      user_ModelObj = user_ModelObj,
                                                                      lower = 1e-20,
                                                                      upper = 20,
                                                                      method = "Brent",
                                                                      control = list(fnscale = -1, trace = 6, maxit = 1),
                                                                      hessian=T))
        nullOut <- capture.output(
          nullLogLikeList[[names(sampleSubsets)[subset]]][[gene]] <- callGetLLByGene(geneEss = 1, useGene = gene,
                                                                                     useSamples = sampleSubsets[[subset]],
                                                                                     sample_effects = sample_effects,
                                                                                     guide_efficiency = guide_efficiency,
                                                                                     user_DataObj = user_DataObj,
                                                                                     user_ModelObj = user_ModelObj))
        if (gene_idx == 1) {
          write(optimOut, file = paste0(tStamp, '_optim_output.txt'))
          write(nullOut,  file = paste0(tStamp, '_null_output.txt'))
        }
      }
      print("finished optimizing one subset:")
      print(names(sampleSubsets)[subset])
      print(sampleSubsets[[subset]])
    }
    return(list(optObjList, nullLogLikeList))
  }

  # -------------------------------- Main Method -------------------------------------

  # Initialize
  geneList <- unique(user_DataObj$guide2gene_map$gene)
  numGenes <- length(geneList)
  if (is.list(subset_genes)) {
    if (any(c(subset_genes[[1]] < 1,
              subset_genes[[2]] > nrow(user_DataObj$dep_counts),
              length(subset_genes) > 2))) {
      cat('Only ',numGenes, ' genes queried.')
      stop('invalid gene subset specified.')
    }
    subsetG <- do.call(seq, subset_genes)
  } else if (is.na(subset_genes[1])) {
    subsetG <- 1:numGenes
  } else if (is.vector(subset_genes)) {
    subsetG <- subset_genes
  } else {
    stop('invalid gene subset argument.')
  }
  guide_efficiency <- rep(1, nrow(user_DataObj$dep_counts))
  guide_essentiality <- rep(1, nrow(user_DataObj$dep_counts))
  sampleSubsets <- list("all"=1:ncol(user_DataObj$dep_counts))
  startEss <- rep(1, length(geneList))
  gene_essentiality <- startEss
  if (is.vector(user_ModelObj$test_sample_subtype_cols, mode='integer')) {
    sampleSubsets$test <- user_ModelObj$test_sample_subtype_cols
    sampleSubsets$ctrl <- sampleSubsets$all[!sampleSubsets$all %in% sampleSubsets$test]
  }

  # 3 - optimize by-gene params;
  # Determine differential depletion in sample subsets by fitting TWO essential
  # parameters per gene
  optFgObj <- optFg(startEss, sampleSubsets, sample_effects,
                    guide_efficiency, geneList[subsetG])
  optObjList <- optFgObj[[1]]
  nullLogLikeList <- optFgObj[[2]]
  gene_essentiality[subsetG] <- sapply(optObjList[['all']], function(i) i$par)

  # OPtional: iterative optimization of sample-specific tau and gene essentailtiy phi.
  # startEss <- gene_essentiality
  # guide_essentiality <- rep(gene_essentiality,
  # times = table(user_DataObj$guide2gene_map$gene))

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
  print(head(gene_results))

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
      print(c(all_scale, test_scale, ctrl_scale))
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
