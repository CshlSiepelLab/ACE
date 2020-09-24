#' Function getMasterLib
#'
#' Process list of master libraries, which may be single or multiple files,
#' which may have one or more replicates.
#' Currently just averaging over total-count-normalized replicates.
#' return data.table of sgRNA x masterlibraries, with column names corresponding
#' to those referenced in sample_masterlib.
#' @import data.table
#' @import plyr
#' @param countList List of data.tables with counts of sequenced masterlibraries.
#' One masterlibrary per data.table; all replicates within one data.table.
#'
#' @return Data.table with sgrna x masterlibs of log frequencies.
#' @export
#'
getMasterLib <- function(countList) {
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
    qmax <- max(masterMean)
    set(mf_avg_dt,
        j = ml_name,
        value = masterMean - qmax - log(sum(exp(masterMean-qmax))))
    print('averaged masterlib replicates for:')
    print(ml_name)
    print(head(mf_avg_dt))
    if (exists('master_freq')) master_freq <- merge(master_freq,
                                                    mf_avg_dt,
                                                    by='sgrna')
    else master_freq <- mf_avg_dt
  }
  setkey(master_freq, sgrna)
  return(master_freq)
}
