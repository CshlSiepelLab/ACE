#' Function writeResObj
#' write results objects
#' @import data.table
#' @param user_ResObj ResObj or list; output from optimizeModelParameters.
#' @return Returns index of written results version.
#' @export
writeResObj <- function(user_ResObj) {
  objName <- deparse(substitute(user_ResObj))
  print("writing results from")
  print(objName)
  v<- 0
  while (paste0(objName, "_", v, '_gene_results.txt') %in% dir()) {
    v <- v + 1
  }
  save(objName, file = paste0(objName, '_v', v))
  write.table(user_ResObj$gene_results,
              file = paste0(objName, "_", v, '_gene_results.txt'),
              row.names = F, quote = F)
  write.table(user_ResObj$diff_genes,
              file=paste0(objName, "_", v, '_diff_genes.txt'),
              row.names = F, quote = F)
  write(paste0(timestamp(),'\nREADME\n',user_ResObj$readMe),
        file = paste0(objName, "_", v, '_readMe.txt'))
  print(paste0(objName, ' written'))
  return(list('gene_results'=paste0(objName, "_", v, '_gene_results.txt'),
              'diff_genes' =paste0(objName, "_", v, '_diff_genes.txt'),
              'readMe' = paste0(objName, "_", v, '_readMe.txt')))
}
