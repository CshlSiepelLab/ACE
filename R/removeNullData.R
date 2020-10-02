#' Function removeNullData
#' 
#' Function to strip out sgRNA & samples with no count data.
#' Called by default within DataObjClass.  Do not specify test_tumor_subtype_cols
#' in DataObjClass
#' until AFTER erroneous samples are removed, or indicing will be messed up.
#' @param user_DataObj Class DataObj, to be stripped. Must be modifiable.
#' @export
removeNullData <- function(user_DataObj) {
  
  # Set local definitions to prevent R check note due to data.table syntax.
  masterlib <- NULL
  sgrna <- NULL
  
  if (is.data.table(user_DataObj$init_counts)) {
    blankInitRows <- which(rowSums(user_DataObj$init_counts)==0)
    blankInitSamples <- which(colSums(user_DataObj$init_counts)==0)
    # update samples in init_counts only.
    if (length(blankInitSamples) >0) {
      warning('Some initial samples have no count data, removing.')
      removeSamples <- blankInitSamples
      user_DataObj$init_counts <- user_DataObj$init_counts[, .SD, 
                                                           .SDcols = -removeSamples]
      # init_total_reads
      user_DataObj$init_total_reads <- user_DataObj$init_total_reads[-removeSamples]
    } else {
      removeSamples <- unique(which(user_DataObj$init_total_reads == 0),
                              which(user_DataObj$cells_infected ==0))  
      if (length(removeSamples) == 0) removeSamples <- F
    }
    # update guides in init_counts only.
    if (length(blankInitRows) > 0) {
      warning('Some guides have no counts in any initial samples, removing.')
      removeGuides <- blankInitRows
      # update init_counts
      user_DataObj$init_counts <- user_DataObj$init_counts[-removeGuides, .SD]
    } else removeGuides <- F
    
  # For experiments with no initial counts. 
  } else {
    # get guides missing from masterlib AND corresponding dep counts.
    sapply(seq_along(user_DataObj$master_counts), function(i) {
      mlibName <- names(user_DataObj$master_counts)[i]
      mlib <- user_DataObj$master_counts[[mlibName]]
      # print(head(mlib))
      samples <- user_DataObj$sample_masterlib[masterlib %in% mlibName, sample]
      blankDepRows <- user_DataObj$dep_counts[, which(rowSums(.SD)==0), 
                                              .SDcols = samples]
      blankDepGuides <- user_DataObj$guide2gene_map$sgrna[blankDepRows]
      blankMasterRows <- mlib[, rowSums(.SD) == 0,
                              .SDcols = -'sgrna']
      blankMasterGuides <- mlib[blankMasterRows, sgrna]
      blankGuides <- blankMasterGuides[blankMasterGuides %in% blankDepGuides]
      # Replace missing masterlib-sample count pairs with NA.
      sapply(blankGuides, function(g) {
        set(user_DataObj$dep_counts,
            i = which(user_DataObj$guide2gene_map$sgrna == g),
            j = samples, value = NA)
      })
    })
    blankGuides <- user_DataObj$dep_counts[, which(rowSums(is.na(.SD))!=0)] 
    
    # remove NA genes and guides.  ACE is for coding regions only.
    naGenes <- which(is.na(user_DataObj$guide2gene_map$gene))
    
    if (length(blankGuides)>0)  removeGuides <- unique(c(blankGuides, naGenes))
    else removeGuides <- F
    
    # remove and flag samples with master libraries with no counts in ANY replicate,
    # and flag master library samples with replicates with no counts.
    lapply(seq_along(user_DataObj$master_counts), function(i) {
      blankS <- user_DataObj$master_counts[[i]][, colSums(.SD)==0, .SDcols = -'sgrna']
      if (all(blankS)) {
        stop(paste0('Master library file',
                    names(user_DataObj$master_counts)[i], 'has no counts.'))
      } else if (any(blankS)) {
        warning('Some master library replicates are blank; removing.')
        print(head(user_DataObj$master_counts[[i]]))
        user_DataObj$master_counts[[i]][, (which(blankS)) := NULL]
        print(head(user_DataObj$master_counts[[i]]))
      }
    })
    if (any(user_DataObj$cells_infected == 0)) {
      removeSamples <- which(user_DataObj$cells_infected == 0)
    } else removeSamples <- F
  }
  
  # Replace missing/dropped guides and samples.
  # test_tumor_subtype_cols - must not be specified in DataObj until AFTER samples
  # removed.
  # user_DataObj$user_DataObj$master_counts - already edited line ~50
  # removeGuides & removeSamples are indices corresponding to count data & 
  # guide2gene_map.
  if (!isFALSE(removeGuides)) {
    warning('Some guides have no counts in any depleted samples or masterlib, removing.')
    tStamp <- paste(unlist(str_split(Sys.time(), ' ')), collapse='_')
    write.table(removeGuides, file = file.path('data',paste0(tStamp,'no_count_guides.txt')))
    if (length(removeGuides) == nrow(user_DataObj$dep_counts)) {
      stop('No valid data submitted.')
    }
    # dep_counts - add pseudocount.
    user_DataObj$dep_counts <- user_DataObj$dep_counts[-removeGuides,.SD, 
                                                       .SDcols = -removeSamples]
    user_DataObj$guide2gene_map <- user_DataObj$guide2gene_map[-removeGuides, ]
  }
  if (!isFALSE(removeSamples)) {
    warning('Some samples have no counts in the masterlib or init samples; removing.')
    write.table(names(user_DataObj$dep_counts)[removeSamples])
    # cells_infected
    user_DataObj$cells_infected <- user_DataObj$cells_infected[-removeSamples]
    # dep_total_reads
    user_DataObj$dep_total_reads <- user_DataObj$dep_total_reads[-removeSamples]
    if (length(removeSamples) == ncol(user_DataObj$dep_counts)) {
      stop('No valid data submitted.')
    }
  }
  
  # Warn if all guides in a depleted sample are 0; could happen if all functional.
  blankDepSamples <- which(colSums(user_DataObj$dep_counts)==0)
  if (length(blankDepSamples) > 0 ) {
    warning('Some depleted samples have no counts in any sgRNA, check input data.')
  }
  # All guide_covars, sample_masterlib info kept; extra info ok.
  
  # # add pseudocount of 0.5 to dep_counts, init_counts, total_counts.
  # user_DataObj$init_counts <- user_DataObj$init_counts + 0.5
  # user_DataObj$dep_counts <- user_DataObj$dep_counts + 0.5
  # user_DataObj$init_total_reads <- user_DataObj$init_total_reads + 
  #   0.5*nrow(user_DataObj$dep_counts)
  # user_DataObj$dep_total_reads <- user_DataObj$dep_total_reads + 
  #   0.5*nrow(user_DataObj$dep_counts)
  # 
  # Return DataObj to allow chaining.
  return(user_DataObj)
}