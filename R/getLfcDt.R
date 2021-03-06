#' getLfcDt: Function to retrieve average log fold change of count data.
#'
#' Takes geometric average across all guides and samples applying to a gene.
#' If negative controls provided, normalizes final averages to mean of negative
#' controls.  Adds a pseudocount of 0.5 to the count data after scaling all
#' counts to 1e7 reads per sample. Returns log2 of dep/init ratio.
#' @param user_DataObj DataObj to analyze.
#' @param isSim Boolean whether data counts simulated; if so, true phi in name.
#' @param use_base_counts Depleted-count-matched data.table to use as the
#' denominator for fold change.
#' @param use_samples Default NA, return results for only a sample subset.
#' @param getCI Default true, return confidence intervals per phi.  Sim only.
#' @param lfcFileName Optional name of output file of lfc values.
#' @param writeFile Boolean default false; should output file be written.
#' @param subset_base Boolean default true; should use_base_counts be subset by
#' use_samples?
#' @param write_log Function to output log file.
#' @param neg_ctrl_file Name of file containing gene names to use as negative
#' controls.  Optional; default is file used by user_DataObj.
#' @export
getLfcDt <- function(user_DataObj, 
                     isSim,
                     use_base_counts=NA, 
                     use_samples=NA, 
                     getCI=T,
                     lfcFileName='lfc', 
                     writeFile = F,
                     subset_base=T,
                     write_log,
                     neg_ctrl_file = NA) {
  # Set local definitions to prevent R check note due to data.table syntax.
  masterlib <- NULL
  gene <- NULL
  raw_avg_lfc <- NULL
  true_gene_param <- NULL
  score <- NULL
  
  genes <- user_DataObj$guide2gene_map$gene
  if (is.na(neg_ctrl_file)) {
    neg_ctrl_file <- user_DataObj$neg_control_file # may still be NA.
  }
  if (any(is.na(use_samples))) use_samples <- 1:ncol(user_DataObj$dep_counts)
  if (!all(use_samples %in% 1:ncol(user_DataObj$dep_counts))) {
    stop('Invalid sample index in use_samples argument.')
  }
  if (any(is.na(use_base_counts))) {
    if (is.data.table(user_DataObj$init_counts)) {
      message('Using initial counts in calculation of LFC.')
      write_log('Using initial counts in calculation of LFC.')
      use_base_counts <- user_DataObj$init_counts
      use_samples_init <- use_samples
    } else {
      stop('provide master library or init counts as base_counts argument')
    }
  } else {
    # check number of genes match.
    if (nrow(user_DataObj$dep_counts) != nrow(use_base_counts) |
        nrow(use_base_counts) != nrow(user_DataObj$guide2gene_map)) {
      stop('use_base_counts has different guides than dep_counts')
    }
    # check number of samples match given subsetting arguments.
    # Case 1: use_samples and use_base_counts match up, and subset_base==F
    # (no column reordering of reference to match use_samples)
    # Case 2: subset_base==T, and use_samples works for use_base_counts.
    # Case 3: use_base_counts is a single column, and needs to be duplicated.
    # Case 4: error.
    if (!subset_base & length(use_samples)==ncol(use_base_counts)) {
      use_base_counts <- use_base_counts
      use_samples_init <- 1:ncol(use_base_counts)
    } else if (subset_base & all(use_samples %in% 1:ncol(use_base_counts))) {
     use_base_counts <- use_base_counts
     use_samples_init <- use_samples
    } else if (ncol(use_base_counts)==1) {
      message('One reference sample; matching to depleted counts.')
      use_base_counts <- use_base_counts[, rep(1, length(use_samples)),
                                         with=F]
      use_samples_init <- 1:ncol(use_base_counts)
    } else {
      message('Dimensions of use_base_counts:', dim(use_base_counts))
      stop('Incorrect number of columns in use_base_counts arg.')     
    }
  }
  init_counts <- use_base_counts[, (use_samples_init), with=F]
  dep_counts <- user_DataObj$dep_counts[, (use_samples), with=F]


  # normalize reads scaling read count total to 10M reads.
  dep_scale <- 1e7/dep_counts[ ,lapply(.SD+.5, sum)]
  init_scale <- 1e7/init_counts[, lapply(.SD+.5, sum)]
  norm_ratio <- dep_scale/init_scale

  lfc_dt <- as.data.table(lapply(seq_along(init_counts), function(i)
    log2((dep_counts[[i]] + 0.5)/(init_counts[[i]]+0.5) * norm_ratio[[i]])))
  names(lfc_dt) <- names(dep_counts)
  avg_lfc <- data.table("gene" = genes,
                        "raw_avg_lfc" = rowMeans(lfc_dt))
  if (!is.na(neg_ctrl_file)) {
    neg_ctrl_genes <- fread(neg_ctrl_file,
                            header = F,
                            sep='')[[1]]
    if (sum(avg_lfc$gene %in% neg_ctrl_genes) > 0) {
      neg_ctrl_lfc <- avg_lfc[gene %in% neg_ctrl_genes, mean(raw_avg_lfc)]
    } else {
      neg_ctrl_lfc <- 0
    }
  } else {
    neg_ctrl_lfc <- 0
  }
  avg_lfc[, 'avg_lfc' := raw_avg_lfc-neg_ctrl_lfc]

  lfc_dt$sgrna <- user_DataObj$guide2gene_map$sgrna
  lfc_dt$gene <- user_DataObj$guide2gene_map$gene
  setcolorder(lfc_dt, neworder = c('sgrna', 'gene'))

  # Write lfc file out.
  if (writeFile) {
    message('lfc_guide_counts saved to file:')
    message(lfcFileName)
    write_log(c('lfc guide counts saved to file:\t', lfcFileName))
    write_log(head(lfc_dt))
    v <- getFileIdx(lfcFileName)
    fileName <- paste0(lfcFileName, v, '.txt')
    write.table(lfc_dt, file = fileName, row.names = F, quote = F, sep = '\t',
                append = F)
  }

  # Average lfc by gene (geometric b/c in log2 space).
  gene_avg_lfc <- avg_lfc[, mean(avg_lfc), by=gene]
  setnames(gene_avg_lfc, 'V1', 'score')
  gene_avg_lfc[, 'method' := 'average_lfc']
  if (isSim) {
    gene_avg_lfc[, "true_gene_param" := tstrsplit(gene, "_", keep = 2)]
    gene_avg_lfc[, true_gene_param := 1- as.numeric(true_gene_param)]

    if (getCI) {
      gene_avg_lfc[, "CI_lower" := sort(score)[round(.2*.N)],by=true_gene_param]
      gene_avg_lfc[, "CI_upper" := sort(score)[round(.9 *.N)], by=true_gene_param]
    }
  }
  return(list('gene_avg_lfc' = gene_avg_lfc,
              'lfc_dt' = lfc_dt))
}
