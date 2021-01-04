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
    cat("arg list lengths and head:\n", file = paste0(tStamp, '_argList.txt'),
        append = F)
    cat('useGene:\n', useGene, '\nuseSamples:\n', useSamples,'\n',
        file = paste0(tStamp, '_argList.txt'), append=T)
    for (i in seq_along(argList)) {
      write(c(argNames[i], 'length: ', length(argList[[i]]), 'head:'),
            append=T,
            file = paste0(tStamp, '_argList.txt'))
      write(argList[[i]][1:10], file = paste0(tStamp, '_argList.txt'),
            append = T)
    }
    save(argList, file = paste0('argList_', tStamp,'.RData'))
  }
  if (isProblem) {
    write_log("NA's in argList, see", tStamp,"_argList.txt")
    stop(cat("NA's in argList, see", tStamp,"_argList.txt"))
  }
}