#' Function writeResObj
#' write results objects
#' @import data.table
#' @param user_ResObj ResObj or list; output from optimizeModelParameters.
#' @return Returns index of written results version.
#' @export
writeResObj <- function(user_ResObj) {
  objName <- deparse(substitute(user_ResObj))
  message("writing results from")
  message(objName)
  tStamp <- paste(unlist(str_split(Sys.time(), ' |:')), collapse='_')
  v<- tStamp
  if (!dir.exists('ACE_output_data')) dir.create('ACE_output_data')
  save(objName, file = file.path('ACE_output_data',paste0(objName, '_v', v)))
  write.table(user_ResObj$gene_results,
              file = file.path('ACE_output_data',paste0(objName, "_", v, '_gene_results.csv')),
              row.names = F, quote = F, sep = ',')
  write.table(user_ResObj$diff_genes,
              file= file.path('ACE_output_data',paste0(objName, "_", v, '_diff_genes.csv')),
              row.names = F, quote = F, sep = ',')
  write(paste0(timestamp(),'\nREADME\n',user_ResObj$readMe),
        file = file.path('ACE_output_data',paste0(objName, "_", v, '_readMe.txt')))
  message(paste0(objName, ' written'))
  return(list('gene_results'=file.path('ACE_output_data', paste0(objName, "_", v, '_gene_results.csv')),
              'diff_genes' = file.path('ACE_output_data', paste0(objName, "_", v, '_diff_genes.csv')),
              'readMe' = file.path('ACE_output_data',paste0(objName, "_", v, '_readMe.txt'))))
}
