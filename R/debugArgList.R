#' Function debugArgList
#' 
#' Check for illegal argument values.
#' @param argList List of args to check.
#' @param argNames Names of specific arguments; default all.
#' @param printArgList Boolean option to print argList.
#' 
debugArgList <- function(argList, argNames, printArgList = F){
  isProblem <-F
  if (any(sapply(argList, function(i) any(is.na(unlist(i)))))) {
    isProblem <- T
    printArgList <- T
  }
  if (any(sapply(argList, function(i) any(is.infinite(unlist(i)))))) {
    isProblem <- T
    printArgList <- T
  }
  if (any(sapply(argList, function(i) any(is.null(unlist(i)))))) {
    isProblem <- T
    printArgList <- T
  }
  
  if (printArgList) {
    tStamp <- paste(unlist(str_split(Sys.time(), ' |:')), collapse='_')
    outFile <- file.path('ACE_output_data', paste0(tStamp, '_argList.txt'))
    cat("arg list lengths and head:\n", file = outFile,
        append = F)
    for (i in seq_along(argList)) {
      write(c(argNames[i], 'length: ', length(argList[[i]]), 'head:'),
            append=T,
            file = outFile)
      write(argList[[i]][1:10], file = outFile,
            append = T)
    }
    # save(argList, file = file.path('ACE_output_data',paste0('argList_', tStamp,'.RData')))
  }
  if (isProblem) {
    print("NA's in argList, see", outFile)
    stop(cat("NA's in argList, see", outFile))
  }
}