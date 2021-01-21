#' Function optFg: By-gene essentiality optimization.
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom utils install.packages installed.packages
#' @import parallel
#' @param startEss Starting gene essentiality.
#' @param sampleSubsets Sample groupings into test, ctrl, and 'all' panels.
#' @param sample_effects By-sample parameters.
#' @param guide_efficiency By-guide efficiency parameter.
#' @param geneList Genes to evaluate sequentially in this function call.
#' @param write_log Function to write to central log file.
#' @param user_DataObj DataObject containing count data to analyze.
#' @param user_ModelObj ModelObject with pre-processed model parameters.
#' @param ncpus Number of threads to use in parallel optimization; if not 
#' specified, will be auto-detected.
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
                  user_ModelObj,
                  ncpus = NA) {
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
  
  optLogLikeList <- list()
  nullLogLikeList <- list()
  ctr <- 0
  progressBar <- txtProgressBar(min = ctr,
                                max = length(geneList),
                                initial=0)
  env <- environment()

  # Optimization Function.
  optFg_OneThread <- function(gene_idx, subset_idx) {
    gene <- geneList[[gene_idx]]
    optimOut <- capture.output(
      optLL <- optim(
        par = startEss[[gene]],
        fn = callGetLLByGene,
        useGene = gene,
        useSamples = sampleSubsets[[subset_idx]],
        sample_effects = sample_effects,
        guide_efficiency = guide_efficiency,
        user_DataObj = user_DataObj,
        user_ModelObj = user_ModelObj,
        # write_log = write_log,
        lower = 1e-20,
        upper = 20,
        method = "Brent",
        control = list(fnscale = -1, trace = 6, maxit = 1),
        hessian=T))
    nullOut <- capture.output(
      nullLL <- callGetLLByGene(
        geneEss = 1, useGene = gene,
        useSamples = sampleSubsets[[subset_idx]],
        sample_effects = sample_effects,
        guide_efficiency = guide_efficiency,
        user_DataObj = user_DataObj,
        user_ModelObj = user_ModelObj))
    return(list('gene' = gene,
                'nullLL' = nullLL,
                'optLL' = optLL,                
                'logMessages' = list('optimOut'=optimOut, 'nullOut'=nullOut)))
  }
  
  # Optimization wrapper for progress reporting.
  progWrapper <- function (gene_idx, subset_idx) {
    currentCount <- get('ctr', envir = env) + 1
    assign('ctr', currentCount, envir=env)
    setTxtProgressBar(get('progressBar', envir=env), currentCount)
    optObjList <- optFg_OneThread(gene_idx, subset_idx)
    return(optObjList)
  }
  
  # using 'parallel' R package to optimize gene parameters in parallel.
  # use half of locally available threads - may conflict with cluster manager
  # assignments unless entire node reserved, or running on a local machine.
  if (!is.numeric(ncpus)) {
    numThreads <- min(floor(detectCores()/2), floor(length(geneList)/10) + 1)
  } else {
    numThreads <- round(ncpus)
  }
  cl <- makeCluster(numThreads)
  # did not help with i386 problem.
  # currentPaths <- .libPaths()
  # clusterExport(cl, varlist = list('currentPaths'), envir = env)
  # clusterEvalQ(cl, .libPaths(currentPaths))
  clusterEvalQ(cl, 
               lapply(c('data.table', 'stats', 'utils', 'Rcpp', 'ACER'), 
                      function(l) {
                 if (! l %in% installed.packages()) install.packages(l)
                 require(l, character.only = T)
               }))
  clusterExport(cl, varlist = list('callGetLLByGene', 
                                   'getLL',
                                   'getLLNoInit',
                                   'geneList',
                                   'sampleSubsets',
                                   'sample_effects',
                                   'startEss',
                                   'guide_efficiency',
                                   'user_DataObj',
                                   'user_ModelObj'),
                envir = env)
  on.exit(stopCluster(cl))
  write_log("Using parallel cluster:")
  write_log(as.character(cl))
  for (subset_idx in seq_along(sampleSubsets)) {
    optList <- parLapply(cl, seq_along(geneList), function(g) {
      progWrapper(g, subset_idx)
    })
    for (optObj in optList) {
      nullLogLikeList[[names(sampleSubsets)[subset_idx]]][[optObj$gene]] <- 
        optObj$nullLL
      optLogLikeList[[names(sampleSubsets)[subset_idx]]][[optObj$gene]] <-
        optObj$optLL
      write_log(c('Optimized gene:', optObj$gene))
      write_log(optObj$logMessages$optimOut)
      write_log(optObj$logMessages$nullOut)
    }
    write_log("finished optimizing gene parameters over one sample subset:")
    write_log(names(sampleSubsets)[subset_idx])
    write_log(sampleSubsets[[subset_idx]])
  }
  # setDefaultCluster(NULL)
  
  gene_effects <- lapply(optLogLikeList[['all']], function(i) i$par)
  names(gene_effects) <- names(optLogLikeList[['all']])
  return(list(optLogLikeList, nullLogLikeList,
              'sample_effects'=sample_effects,
              'guide_efficiency' = guide_efficiency,
              'gene_effects' = gene_effects))
}