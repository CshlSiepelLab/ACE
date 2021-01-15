#' Function optFg: By-gene essentiality optimization.
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @param startEss Starting gene essentiality.
#' @param sampleSubsets Sample groupings into test, ctrl, and 'all' panels.
#' @param sample_effects By-sample parameters.
#' @param guide_efficiency By-guide efficiency parameter.
#' @param geneList Genes to evaluate sequentially in this function call.
#' @param write_log Function to write to central log file.
#' @param user_DataObj DataObject containing count data to analyze.
#' @param user_ModelObj ModelObject with pre-processed model parameters.
# assume each gene is independent, and optimize over a single gene at a time.
# while this single parameter optimization is very fast, this for loop *can* be parallelized.
# produces list of optim objects, one per sample subgroup.
# Assumes phi's positive selection impact capped at 20-fold (upper optim bound).
#mu_1 = lambda1[s]*n
#mu_2 = lambda2[s]* mu_1 * (1 - eps + eps * phi)
# (see src/getLL.cpp for the actual math).
optFg <- function(startEss, sampleSubsets, sample_effects,
                  guide_efficiency, geneList, write_log,
                  user_DataObj,
                  user_ModelObj) {
  # Check for empty input.
  if (length(sampleSubsets)==0 | any(is.na(sampleSubsets))) {
    write_log('Error: Sample subset submitted to optFg is empty or NA')
    write_log(sampleSubsets)
    stop('Sample subset submitted to optFg is empty or NA')
  }
  if (length(geneList)==0 | any(is.na(geneList))) {
    write_log('Error: Gene list submitted to optFg is empty or NA')
    write_log('---')
    write_log(head(geneList))
    write_log('---')
    stop('Error: Gene list submitted to optFg is empty or NA')
  }
  
  optObjList <- list()
  nullLogLikeList <- list()
  progressBar <- txtProgressBar(min = 0,
                                max = length(sampleSubsets) * length(geneList),
                                initial=0)
  
  for (subset in seq_along(sampleSubsets)) {
    for (gene_idx in seq_along(geneList)) {
      gene <- geneList[[gene_idx]]
      write_log(c("optimizing parameters for gene: ", gene, '\n'))
      optimOut <- capture.output(
        optObjList[[names(sampleSubsets)[subset]]][[gene]] <- optim(par = startEss[[gene]],
                                                                    fn = callGetLLByGene,
                                                                    useGene = gene,
                                                                    useSamples = sampleSubsets[[subset]],
                                                                    sample_effects = sample_effects,
                                                                    guide_efficiency = guide_efficiency,
                                                                    user_DataObj = user_DataObj,
                                                                    user_ModelObj = user_ModelObj,
                                                                    write_log = write_log,
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
                                                                                   user_ModelObj = user_ModelObj,
                                                                                   write_log = write_log))
    }
    write_log("finished optimizing one subset:")
    write_log(names(sampleSubsets)[subset])
    write_log(sampleSubsets[[subset]])
    setTxtProgressBar(progressBar, gene_idx + gene_idx * (subset - 1))
  }
  gene_effects <- lapply(optObjList[['all']], function(i) i$par)
  names(gene_effects) <- geneList
  return(list(optObjList, nullLogLikeList,
              'sample_effects'=sample_effects,
              'guide_efficiency' = guide_efficiency,
              'gene_effects' = gene_effects))
}