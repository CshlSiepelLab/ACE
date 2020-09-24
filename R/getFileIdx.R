#' Function getFileIdx
#' 
#' Returns run index.
#' @param baseFileName Vector of names of file to check for, for earlier runs 
#' with conflicting names.  Appended with '#.txt'
#' @examples 
#' getFileIdx('countData_sim_10_genes_phi_G1_v')
#' @export
getFileIdx <- function(baseFileName) {
  fileIdx <- 0
  print('checking current directory for previous results:')
  print(getwd())
  fileNames <- dir()
  fullFile <- paste0(baseFileName, fileIdx, '.txt')
  while (any(fullFile %in% fileNames)) {
    fileIdx <- fileIdx + 1
    fullFile <- paste0(baseFileName, fileIdx, '.txt')
    print(fullFile)
  }
  return(fileIdx)
}
    
  