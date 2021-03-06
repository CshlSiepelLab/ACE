#' Function deriveSampleScaling
#' create init_scaling and dep_scaling. 
#' Misspecification will cause our 'neutral' model to
#' actually be depleted, reducing all essentiality scores.
#' The expected scaling factor is the ratio of read counts to their
#' corresponding prior, with the median taken across sgrna in a sample.  
#' This method assumes the median sgrna is nonessential.
#' This is also approximated by the ratio of the total sample read counts
#' to the total expected read counts.
#' If indicated, the scaling factor will be calculated only using negative
#' control sgrna (recommended).
#' If the master library is used as the prior, the expected read count is
#' the product of the master library sgrna frequency and the total number of
#' infected cells.
#' If the master library sgrna frequency is not used, the expected frequency for
#'  each guide is determined by the 
#' average percent abundance for that guide across all (mean-normalized) initial
#' samples. 
#' @param user_DataObj all data; a DataObj.
#' @param master_freq_dt normalized masterlibrary frequencies (log)
#' @param use_master_library Boolean
#' @param use_neg_ctrl Boolean
#' @param neg_ctrls Vector of gene name strings.
#' @param write_log Function to write messages to log.
#' @importFrom stats sd median
deriveSampleScaling <- function(user_DataObj,
                                master_freq_dt,
                                use_master_library,
                                use_neg_ctrl,
                                neg_ctrls,
                                write_log) {
  # locally binding variable names used for DataObj columns.
  sgrna <- NULL
  # sanity check.
  if (!any(user_DataObj$dep_counts[, colSums(.SD)] > 0)) {
    stop('No counts in any depleted samples, provide data for negative controls.')
  }
  
  if (use_neg_ctrl) {
    useGuides <- which(user_DataObj$guide2gene_map$gene %in% neg_ctrls)
  } else {
    useGuides <- 1:nrow(user_DataObj$dep_counts)
  }

  # log of expected initial guide abundance.
  base_counts <- list()
  if (use_master_library & is.data.table(user_DataObj$sample_masterlib)) {
    write_log('Scaling samples relative to matched masterlibrary')
    message('Scaling samples relative to matched masterlibrary')
    for (i in seq_along(names(user_DataObj$dep_counts))) {
      sampleName <- names(user_DataObj$dep_counts)[i]
      useMasterlib <- unlist(user_DataObj$sample_masterlib[sample==sampleName,
                                                           masterlib])
      masterlib <- master_freq_dt[user_DataObj$guide2gene_map$sgrna,
                                  ][[useMasterlib]]
      base_counts[[i]] <- masterlib[useGuides]
    }
    # If no neg ctrls being used, check for over-depleted libraries.
    if (!use_neg_ctrl) {
      matchedMlibDT <- as.data.table(lapply(names(user_DataObj$dep_counts),
                                            function(i) {
                                              useMasterLib <- user_DataObj$sample_masterlib[sample == i,
                                                                                            masterlib]
                                              master_freq_dt[sgrna %in% user_DataObj$guide2gene_map$sgrna,
                                                             (useMasterLib), with=F]
                                            }))
      numGuides <- nrow(user_DataObj$dep_counts)
      deplete_50p_counts <- rbind(matchedMlibDT[sample.int(numGuides,
                                                           size=.5*numGuides,
                                                           replace=F), .SD],
                                  matrix(0, nrow = .5*numGuides,
                                         ncol = ncol(user_DataObj$dep_counts)),
                                  use.names = F)
    }
  } else if (is.data.table(user_DataObj$init_counts)) {
    write_log('Scaling samples relative to mean abundance in initial samples.')
    message('Scaling samples relative to mean abundance in initial samples.')
    # Assuming no counts are depleted in the initial population.
    for (j in seq_along(names(user_DataObj$init_counts))) {
      guideMeans <- user_DataObj$init_counts[useGuides, exp(rowMeans(log(.SD+0.5)))]
      normCounts <- user_DataObj$init_counts[useGuides,(0.5+.SD)/sapply(seq_along(.SD),function(i)
                                               median((0.5+user_DataObj$init_counts[[i]])/guideMeans[i]))]
      # get average percent abundance of each guide across all samples: guideMeans/sum(guideMeans)
      # assuming a common masterlibrary for all.
      base_counts[[j]] <- log(guideMeans/sum(guideMeans))
      # get by-sample percent abundance in initial counts, based on normalized scores:
      # get by-guide average ratio of normalized scores to total counts.
      # base_counts[[j]] <- rowMeans(normCounts[useGuides, log(.SD)] -
                                     # user_DataObj$init_counts[,log(colSums(0.5+.SD))])
    }
    # If no neg ctrls being used, check for over-depleted libraries.
    if (!use_neg_ctrl) {
      numGuides <- nrow(user_DataObj$init_counts)
      deplete_50p_counts <- rbind(user_DataObj$init_counts[
        sample.int(numGuides, size=.5*numGuides, replace=F), .SD],
        matrix(0, nrow = .5*numGuides,
               ncol = ncol(user_DataObj$init_counts)),
        use.names = F)
    }
  } else {
    stop('Must use the master library or initial abundances to scale samples.')
  }
  
  #check for overdispersion if no neg ctrls used.
  # compare sd of depleted counts to the most variable
  # masterlibrary counts, artificially depleted (50%)
  if (!use_neg_ctrl) {
    ref_scaling_ratio <- unlist(deplete_50p_counts[, lapply(.SD, function(i)
      sd(i/mean(i)))])
    dep_scaling_ratio <- unlist(user_DataObj$dep_counts[,lapply(.SD, function(i)
      sd(i/mean(i)))])
    if(any(dep_scaling_ratio > ref_scaling_ratio)) {
      message('Excess of guides seem depleted; recommend using negative controls.')
      message('Using controls?: ', use_neg_ctrl, "\n")
      write_log('Excess of guides seem depleted; recommend using negative controls.')
      write_log(c('Using controls?: ', use_neg_ctrl, "\n"))
    }
  }  
  
  # Calculate scaling relative to base counts.
  if (is.data.table(user_DataObj$init_counts)) {
    init_scaling <- sapply(seq_along(user_DataObj$init_counts),
                               function(i) {
                                 exp(median(log(0.5 + user_DataObj$init_counts[[i]][useGuides])-
                                              base_counts[[i]]) -
                                       log(user_DataObj$cells_infected[i]))})
    if (any(sapply(init_scaling, is.na))) {
      naIdx <- which(is.na(init_scaling))
      write_log("ERROR: NA's produced in calculation of init_scaling")
      write_log(c('using guide indices: ',useGuides[1:10]))
      if (use_neg_ctrl) write_log(c('10 neg ctrls: ',neg_ctrls[1:10]))
      write_log(c('init samples at na: ',names(user_DataObj$init_counts)[naIdx]))
      write_log(c('init counts in na samples: ',
                  user_DataObj$init_counts[useGuides[1:10],naIdx, with=F]))
      write_log(c('cells infected in na samples: '))
      write_log(user_DataObj$cells_infected[naIdx])
      write_log(c('base counts: ',
                  sapply(base_counts[naIdx], function(i) i[1:10])))
      testExp <- names(user_DataObj$dep_counts)[[1]]
      write_log('Sample base count calculation with:')
      write_log(testExp)
      write_log(c('useMasterLib: ', unlist(user_DataObj$sample_masterlib[sample==testExp,
                                                                         masterlib])))
      write_log(c('masterlib ', head(master_freq_dt[user_DataObj$guide2gene_map$sgrna,][[
        useMasterlib]][useGuides])))
      stop('NA in init_scaling.')
    }
  } else init_scaling <- NA
  
  dep_scaling <- sapply(seq_along(user_DataObj$dep_counts),
                        function(j) {
                          exp(median(log(0.5 + user_DataObj$dep_counts[[j]][useGuides])-
                                       base_counts[[j]]) -
                                log(user_DataObj$cells_infected[j]))})
  
  if (any(sapply(dep_scaling, is.na))) {
    write_log('ERROR: NA in scaling parameter.')
    naIdx <- which(is.na(dep_scaling))
    write_log(names(user_DataObj$dep_counts)[naIdx])
    write_log(master_freq_dt[1:10, .SD])
    write_log(user_DataObj$dep_counts[useGuides[1:10], (naIdx), with=F])
    stop('NA in dep_scaling')
  }


  return(list('init_scaling'=init_scaling, 'dep_scaling'=dep_scaling))
}
