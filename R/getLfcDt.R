#' getLfcDt: Function to retrieve average log fold change of count data.
#'
#' Takes geometric average across all guides and samples applying to a gene.
#' If negative controls provided, normalizes final averages to mean of negative
#' controls.  Adds a pseudocount of 0.5 to the count data after scaling all
#' counts to 1e7 reads per sample. Returns log2 of dep/init ratio.
#' @param user_DataObj DataObj to analyze.
#' @param isSim Boolean whether data counts simulated; if so, true phi in name.
#' @param master_freq_dt Data.table of master library abundances if initial read counts absent.
#' @param use_base_counts Depleted count -matched data.table to use as the
#' denominator for fold change.
#' @param use_samples Default NA, return results for only a sample subset.
#' @param getCI Default true, return confidence intervals per phi.  Sim only.
#' @param lfcFileName Optional name of output file of lfc values.
#' @param writeFile Boolean default false; should output file be written.
#' @export
getLfcDt <- function(user_DataObj, 
                     isSim,
                     master_freq_dt = NA,
                     use_base_counts=NA, 
                     use_samples=NA, 
                     getCI=T,
                     lfcFileName='lfc', 
                     writeFile = F) {
  # Set local definitions to prevent R check note due to data.table syntax.
  masterlib <- NULL
  gene <- NULL
  raw_avg_lfc <- NULL
  true_gene_param <- NULL
  score <- NULL
  
  genes <- user_DataObj$guide2gene_map$gene
  neg_ctrl_file <- user_DataObj$neg_control_file
  if (is.na(use_samples[1])) use_samples <- 1:ncol(user_DataObj$dep_counts)
  if (is.na(use_base_counts)) {
    if (is.data.table(user_DataObj$init_counts)) {
      print('Using initial counts in calculation of LFC.')
      use_base_counts <- user_DataObj$init_counts
    } else {
      print('Using master library in calculation of LFC.')
      if (!is.data.table(master_freq_dt)) stop('Provide initial or masterlib counts to lfc.')
      # TODO: parse user_DAtaObj$master_counts
      base_counts <- list()
      for (i in seq_along(names(user_DataObj$dep_counts))) {
        sampleName <- names(user_DataObj$dep_counts)[i]
        useMasterlib <- unlist(user_DataObj$sample_masterlib[sample==sampleName,
                                                             masterlib])
        base_counts[[i]] <- master_freq_dt[user_DataObj$guide2gene_map$sgrna,
                                    ][[useMasterlib]]
      }
      print(head(base_counts[[1]]))
      use_base_counts <- as.data.table(base_counts)
      print(head(use_base_counts))
    }
  } else if (ncol(use_base_counts) != ncol(user_DataObj$dep_counts)) {
      if (ncol(use_base_counts)==1) {
        use_base_counts <- as.data.table(rep(use_base_counts[[1]],
                                             ncol(user_DataObj$dep_counts)))
      }
      else stop('incorrect dimensions of submitted base counts.')
  }
  init_counts <- use_base_counts[, (use_samples), with=F]
  dep_counts <- user_DataObj$dep_counts[, (use_samples), with=F]

  # print(head(init_counts))
  # print(dim(init_counts))
  # print(head(dep_counts))
  # print(dim(dep_counts))

  # normalize reads scaling read count total to 10M reads.
  dep_scale <- 1e7/dep_counts[ ,lapply(.SD+.5, sum)]
  init_scale <- 1e7/init_counts[, lapply(.SD+.5, sum)]
  norm_ratio <- dep_scale/init_scale
  # print(norm_ratio)

  lfc_dt <- as.data.table(lapply(seq_along(init_counts), function(i)
    log2((dep_counts[[i]] + 0.5)/(init_counts[[i]]+0.5) * norm_ratio[[i]])))
  names(lfc_dt) <- names(dep_counts)
  avg_lfc <- data.table("gene" = genes,
                        "raw_avg_lfc" = rowMeans(lfc_dt))
  if (!is.na(neg_ctrl_file)) {
    neg_ctrl_genes <- fread(neg_ctrl_file)[[1]]
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
    print('lfc_guide_counts written')
    print(head(lfc_dt))
    v <- getFileIdx(lfcFileName)
    fileName <- paste0(lfcFileName, v, '.txt')
    write.table(lfc_dt, file = fileName, row.names = F, quote = F, sep = '\t',
                append = F)
  }
  # print('averaging lfc file')
  # print(head(avg_lfc))

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
  return(gene_avg_lfc)
}
