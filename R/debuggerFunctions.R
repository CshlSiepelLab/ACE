#' Function debugRes
#' 
#' Returns list of contents of ResObj.  For debugging with testthat.
#' Written April 2019
#' Revised July 2019 to include depletion-only framework & other tweaks.
#' @param user_ResObj Object to analyze.
#' @param user_ModelObj ModelObj used to create ResObj
#' @param user_DataObj DataObj used to create ModelObj & ResObj.
#' @return List of first 10 entries of objects within user_ResObj.
#' @export
debugRes <- function(user_ResObj, user_ModelObj=NA, user_DataObj=NA) {
  childObj <- c('gene_results', 'diff_genes', 'readMe')
  debugList <- list()
  for (c in childObj) {
    debugList[[c]] <- ifelse(is.vector(user_ResObj[[c]]),
                             user_ResObj[[c]][1:10],
                             user_ResObj[[c]][1:10, ])
  }
}

#' Function debugData
#' 
#' Returns list of contents of DataObj.  For debugging with testthat.
#' These are to test for programmer, rather than user, error.
#' @param user_DataObj DataObj used to create ModelObj & ResObj.
#' @return List of first 10 entries, with $flag argument for unusual entries.
#' @export
debugData <- function(user_DataObj, printAll = F) {
  debugList <- list()
  if (is.data.table(user_DataObj$init_counts)) {
    if (any(colMeans(user_DataObj$init_counts)>500)) {
      debugList$flag <- "unusually high mean count for initial reads"
      debugList$init_counts <- user_DataObj$init_counts[1:10, .SD]
    }  
  }
  
  requiredChildObj <- c("cells_infected", "dep_total_reads", "dep_counts")
  childObj <- c('sample_masterlib', 'cells_infected', 'dep_total_reads',
                'init_total_Reads', 'guide2gene_map', 'guide_covars',
                'test_tumor_subtype_cols', 'master_counts', 'dep_counts',
                'init_counts')
  # general test for numeric and NA (some experimental designs purposefully NA).
  isNum <- sapply(requiredChildObj, function(i) all(sapply(user_DataObj[[i]], is.numeric)))
  isNa <- sapply(requiredChildObj, function(i) any(is.na(user_DataObj[[i]])))
  if (any(!isNum)) {
    debugList$flag <- c(debugList$flag, 'non-numeric data entered in mandatory arguments.')
    printAll <- T
  }
  if (any(isNa)) {
    debugList$flag <- c(debugList$flag, "NA's detected in mandatory arguments.")
    printAll <- T
  }
  
  # make sure lengths matching up, have option for no initial counts.
  if (!all.equal(length(user_DataObj$cells_infected),
                 length(user_DataObj$dep_total_reads),
                 ncol(user_DataObj$dep_counts))){
    debugList$flag <- c(debugList$flag, 'vector lengths of mandatory args faulty!')
    printAll <- T
  }
  
  for (c in childObj) {
    if (is.data.table(user_DataObj[[c]]) | is.data.frame(user_DataObj[[c]])) {
      debugList[[c]] <- user_DataObj[[c]][1:10,]
    } else if (is.list(user_DataObj[[c]])) {
      debugList[[c]] <- lapply(user_DataObj[[c]], function(i) i[1:10])
    } else {
      debugList[[c]] <- user_DataObj[[c]][1:10]
    }
  }
  
  if (printAll) {
    print(names(debugList))
    sapply(debugList, print)
  }
}

#' Function debugModel
#' 
#' Returns list of contents of DataObj.  For debugging with testthat.
#' These are to test for programmer, rather than user, error.
#' @param user_DataObj DataObj used to create ModelObj & ResObj.
#' @param user_ModelObj ModelObj from DataObj.
#' @param extraObj If you wanna throw in an extra object, go ahead.
#' @return List of first 10 arg indices, with $flag argument for unusual entries.
#' @export
debugModel <- function(user_ModelObj, user_DataObj, extraObj = c(), printAll = F) {
  debugList <- list()
  if (is.data.table(user_DataObj$init_counts) & 
      any(is.na(user_ModelObj$init_scaling) | length(user_ModelObj$init_scaling)==0)) {
    debugList$flag <- c("init counts present, but init_scaling not calculated.")
    printAll <- T
    debugList$init_counts <- user_DataObj$init_counts[1:10,]
    debugList$init_scaling <- user_ModelObj$init_scaling
  }
  childObj <- c('dep_scaling', 'init_scaling', 'mean_var_model', 
                'unobserved_infected_cell_values', 'master_freq_dt',
                'guide_covar', 'mean_var_model_params', extraObj)

  for (c in childObj) {
    if (is.data.table(user_ModelObj[[c]]) | is.data.frame(user_ModelObj[[c]])) {
      debugList[[c]] <- list(dim(user_ModelObj[[c]]),
                             user_ModelObj[[c]][1:10,])
    } else if (is.list(user_ModelObj[[c]])) {
      debugList[[c]] <- lapply(user_ModelObj[[c]], function(i) {
        list(length(i), i[1:10])})
    } else {
      debugList[[c]] <- list(length(user_ModelObj[[c]]), user_ModelObj[[c]][1:10])
    }
  }
  
  if (printAll) {
    print('printing child object name, dimensions/length, and first 10 indices.')
    for (i in names(debugList)) {
      print(i)
      print(debugList[[i]])
    }
  }
}

#' Function debugCpp
#' 
#' For debugging based on debuggging output file from optimizeModelParameters.
#' @param xDataObj DataObj used to create ModelObj & ResObj.
#' @param xModelObj Created by xDataObj
#' @param cppOutputFile output file.
#' @param numS Number of samples; use to calculate actual number of guides from vector.
#' @return List of first 10 arg indices, with $flag argument for unusual entries.
#' @export
debugCpp <- function(cppOutputFile, xDataObj, xModelObj, numS=NA) {
  # Data.table column names from DataObj & ModelObj
  gene <- NULL
  
  geneOptimOutput <- fread(cppOutputFile)
  if (is.numeric(numS)) {
    numSnumGuides <- numS * xDataObj$guide2gene_map[,.N,by=gene]$N[[1]]
  } else {
  numSnumGuides <- ncol(xDataObj$dep_counts) * xDataObj$guide2gene_map[,.N,by=gene]$N[[1]]
  }
  # data.table col. names.
  maxMean <- NULL
  counts <- NULL
  maxNsg <- NULL
  varType <- NULL
  ll <- NULL
  optimIter <- NULL 
  depInitRatio <- NULL
  Sample <- NULL 
  depScaleN <- NULL
  initScaleN <- NULL
  
  geneOptimOutput[, 'll' := rep(c('inf', 'init', 'dep'), nrow(geneOptimOutput)/3)]
  names(geneOptimOutput) <- c('maxMean', 'counts', 'maxNsg', 'varType')
  geneOptimOutput[,'optimIter' := rep(1:(nrow(geneOptimOutput)/(numSnumGuides*3)), 
                                      each=numSnumGuides*3)]
  geneOptimOutput[, optimIter := as.factor(optimIter)]
  #plotOpt <- ggplot(geneOptimOutput, aes(x=varType, colour = optimIter, y=maxNsg, alpha = 0.2)) + geom_point()
  geneOptimOutput[, 'depInitRatio' := rep(geneOptimOutput[varType=='init', maxNsg]/geneOptimOutput[varType=='dep',
                                                                                                   maxNsg], each=3)]
  #plotDepOpt <- ggplot(geneOptimOutput, aes(x=optimIter, y=depInitRatio)) + geom_point()
  
  print('mean nsg ratio at final iteration is:')
  print(mean(geneOptimOutput[optimIter==nrow(geneOptimOutput)/(numSnumGuides*3),
                             depInitRatio]))
  geneOptimOutput[, 'Sample':= rep(rep(1:ncol(xDataObj$dep_counts), each=3), 
                                   nrow(geneOptimOutput)/(3*ncol(xDataObj$dep_counts)))]
  geneOptimOutput[, 'depScaleN' := xModelObj$dep_scaling[Sample]*maxNsg]
  geneOptimOutput[, 'initScaleN' := xModelObj$init_scaling[Sample]*maxNsg]
  print('mean ratio of maxN * dep_scaling at first and last iteration is:')
  print(mean(geneOptimOutput[varType=='dep' & optimIter==1, counts/depScaleN]))
  print(mean(geneOptimOutput[varType=='dep' & optimIter==nrow(geneOptimOutput)/(numSnumGuides*3),
                             counts/depScaleN]))
  return(geneOptimOutput)
}
