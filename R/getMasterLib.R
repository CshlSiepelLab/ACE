#' Function getMasterLib
#'
#' Process list of master libraries, which may be single or multiple files,
#' which may have one or more replicates.
#' Currently just averaging over total-count-normalized replicates.
#' Returns data.table of sgRNA x masterlibraries, with column names corresponding
#' to those referenced in sample_masterlib.
#' @import data.table
#' @import plyr
#' @param countList List of data.tables with counts of sequenced masterlibraries.
#' One masterlibrary per data.table; all sequencing replicates within one data.table.
#' @param write_log Function to output messages to log file.
#'
#' @return Data.table with sgrna x masterlibs of log frequencies.
#' @export
#'
getMasterLib <- function(countList, write_log) {
  # Set local definitions to prevent R check note due to data.table syntax.
  sgrna <- NULL
  
  # assign each master library to be average of reps.
  for (i in seq_along(countList)) {
    mf_avg_dt <- data.table('sgrna' = countList[[i]]$sgrna)
    ml_name <- names(countList)[i]
    # normalize by mean of master freqs across total-count-normalized replicates.
    # assumes all replicates within a single masterlib count file.
    # pools each count file into a single column of log frequencies.
    masterMean <- countList[[i]][, rowMeans(log(.SD+0.5)), .SD = -c('sgrna')]
    if (any(is.na(masterMean))) write_log(c('NAs in masterMean', masterMean))
    qmax <- max(masterMean)
    set(mf_avg_dt,
        j = ml_name,
        value = masterMean - qmax - log(sum(exp(masterMean-qmax))))
    write_log('averaged masterlib replicates for:')
    write_log(ml_name)
    write_log(head(mf_avg_dt))
    write_log(head(countList))
    write_log('countLIst names are:')
    write_log(names(countList[[i]]))
    if (exists('master_freq')) master_freq <- merge(master_freq,
                                                    mf_avg_dt,
                                                    by='sgrna')
    else master_freq <- mf_avg_dt
  }
  setkey(master_freq, sgrna)
  if (any(is.na(master_freq))) {
    write_log('ERROR: NAs in master library frequencies.')
    stop('NAs in master library frequencies.')
  }
  return(master_freq)
}
