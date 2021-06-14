#' DataObjClass: R6 Environment Object for CRISPR Data
#' 
#' @description 
#' This function is an R6 object class to create an environment for data analysis.
#' Use to load all data associated with a single-target negative selection
#' CRISPR screen; stores all data with minimal processing, but checks for
#' errors and data inconsistencies.
#' 
#' @details 
#' All data objects in this class currently public, permitting operations without
#' copying over the entire dataset in function calls.
#' Normalization results and other preprocessing performed by ModelObjClass.
#'
#' @docType class
#' @import R6
#' @importFrom R6 R6Class
#' @import data.table
#' @return DataObj(R6) of all data used in experiment.  When applicable, raw
#'          count data is returned.  Normalization should happen as ModelObj.
#' @export
DataObj <- R6Class("DataObj",
                   private = list(
                     write_log = 'function'),
                   public = list(
                     #' @field init_counts data.table of initial read counts.
                     init_counts = NA,
                     #' @field dep_counts data.table of depleted read counts.
                     dep_counts = NA,
                     #' @field master_counts data.table of read counts from master library.
                     master_counts = NA,
                     #' @field neg_control_file Name of file containing negative control genes.
                     neg_control_file = NA,
                     #' @field guide_covars File name for sgRNA efficiency prediction features.
                     guide_covars = NA,
                     #' @field guide2gene_map Index for gene targeted by each sgRNA.
                     guide2gene_map = NA,
                     #' @field init_total_reads Total reads per initial sample.
                     init_total_reads = NA,
                     #' @field dep_total_reads Total reads per depleted sample.
                     dep_total_reads = NA,
                     #' @field cells_infected Total number of cells infected per sample.
                     cells_infected = NA,
                     #' @field sample_masterlib Index for which masterlibrary went to which sample.
                     sample_masterlib = NA,
                     #' @field sample_annotations Sample annotations; used for partitioning test/ctrl.
                     sample_annotations = NA,

                     #' @description
                     #' Create a new DataObj
                     #' @param masterFiles Vector of file names of master libraries.  1 column sgrna.
                     #' @param countFile Text file with columns sgrna-name, gene, init_sample_1, dep_sample_1, etc.
                     #' @param negCtrlFile Text file with first column containing names of negative control genes.
                     #' @param cellsInfected Integer of total number of cells infected by masterlibrary in masterFile;
                     #' defaults to 1e3 cells per guide in master library.
                     #' @param guideCovarFile File with features to use in guide efficiency modelling.
                     #' @param sampleInfoFile File with sample initial reads info and sample types (e.g. test/ctrl)
                     #' @param totalReads Optional numeric value or vector with by-sample total read counts.
                     #' @param hasInitSeq Default True; are counts for initial infected guides
                     #'                   present in count file?
                     #' @param sampleMasterInfoFile File mapping samples to masterlibrary files.
                     #' @param useSamples Optional.  Sample indices in countFile to analyze.
                     #' Subsets data stored in object.
                     #' Only depleted sample indices must be given; any initial samples will be auto-matched.
                     initialize = function(masterFiles=NA,
                                           countFile,
                                           negCtrlFile = NA,
                                           cellsInfected = NA,
                                           guideCovarFile = NA,
                                           sampleInfoFile = NA,
                                           totalReads = NA,
                                           hasInitSeq = T,
                                           sampleMasterInfoFile = NA,
                                           useSamples = NA) {
                       
                       # Create log file for info, warnings, and errors.
                       if (!dir.exists('ACE_output_data')) dir.create('ACE_output_data')
                       tStamp <- paste(unlist(str_split(Sys.time(), ' |:')), collapse='_')
                       log_file <- file(file.path('ACE_output_data',
                                                  paste0('ACE_DataObj_log_', tStamp,
                                                         '.txt')),
                                        open='w+')
                       on.exit(close(log_file), add=T)
                       
                       #'  Write messages and warnings to log file.
                       #'  @param message_vector String or table to write.
                       private$write_log <- function(message_vector) {
                         if (is.atomic(message_vector)) {
                           cat(message_vector,'\n', file=log_file, append=T)
                         } else {
                           suppressWarnings(write.table(message_vector, file = log_file, 
                                                        col.names = T, append = T,
                                                        row.names = F))
                         }
                       }
                       
                       # Read optional guide and sample information files
                       if (!is.na(guideCovarFile)) {
                         guide_covars <- fread(guideCovarFile)
                         if (ncol(guide_covars) != 2) {
                           private$write_log("submitted sgrna covariate file is:")
                           private$write_log(head(guide_covars))
                           private$write_log(c("Error:seqFile must have form [sgrna, sequence] with header;",
                                            "may add more covariates later."))
                           stop("seqFile must have form [sgrna, sequence] with header;
                                            may add more covariates later.")
                         }
                       }
                       # have comma-delimited 2nd column, with sample annotations.
                       if (!is.na(sampleInfoFile)) {
                         tryCatch(rawSampleInfo <- fread(sampleInfoFile, header=T),
                                                   warning = function(i) simpleError('Invalid sampleInfoFile'))
                         if (is.null(dim(rawSampleInfo))) stop('Single column sample info file.')
                         if (ncol(rawSampleInfo) == 2) {
                           names(rawSampleInfo) <- c("sample_name", "sample_subtype")
                         } else if (ncol(rawSampleInfo == 3)) {
                           private$write_log('Assuming central column is sample group name.')
                           private$write_log(names(rawSampleInfo))
                           names(rawSampleInfo) <- c('sample_name', 'sample_group', 'sample_subtype')
                         } else {
                           private$write_log("submitted sample metadata file is:")
                           private$write_log(head(rawSampleInfo))
                           private$write_log(c("Error: sample meta data must be ",
                                               "in form [sample_name, sample_annotations]"))
                           stop("sample meta data must be in form [sample_name, sample_annotations]")
                         }
                         self$sample_annotations <- rawSampleInfo
                       }
                       
                       # Load Count Data
                       private$write_log('Reading count file.')
                       dupRawData <- tryCatch(fread(countFile, header = T),
                                              warning = function(i) simpleError('Invalid countData file'),
                                              error = function(i) simpleError('Invalid countData file'))
                       message('Count file read.')
                       # check for duplicates, keep first instance.
                       if (any(duplicated(dupRawData[[1]]))) {
                         private$write_log('using first column as UNIQUE guide labels.')
                         private$write_log(head(dupRawData))
                         private$write_log('Some guides listed twice; discarding both.')
                         message('Some guides listed twice; discarding both.')
                         private$write_log(c('See dup_guides.txt for all duplicated counts;',
                                 'total are:'))
                         private$write_log(sum(duplicated(dupRawData[[1]])))
                         write.table(dupRawData[duplicated(dupRawData[[1]]),],
                                     file = file.path('ACE_output_data','dup_guides.txt'),
                                     quote = F)
                         rawData <- unique(dupRawData, by = 1)
                       } else {
                         rawData <- dupRawData
                       }
                       setorderv(rawData, names(rawData)[2])
                       
                       # check for non-numeric data; NA permitted.
                       if (any(sapply(rawData[,-(1:2), with=F],
                                      function(i) any(!is.na(i) & !is.numeric(i))))) {
                         private$write_log('Non-numeric and non-NA data in raw count file.')
                         invalidCol <- which(sapply(rawData[,-(1:2), with=F],
                                      function(i) any(!is.na(i) & !is.numeric(i)))) 
                         private$write_log(head(rawData[,(invalidCol)+2, with=F]))
                         private$write_log('Often caused be mis-ordered columns.')
                         stop("non-numeric and not-NA data in raw count file. Columns ordered?")
                       }
                       # split into initial & dep counts if applicable.
                       if (hasInitSeq) {
                         if (ncol(rawData) < 4 || ncol(rawData)%%2 == 1) {
                           private$write_log("input count file is:")
                           private$write_log(head(rawData))
                           private$write_log(c('raw counts MUST be in form',
                                               '[sgrna, gene, s1_init, s1_dep,...]',
                                               ' with header.'))
                           stop("raw counts MUST be in form [sgrna, gene, s1_init, s1_dep, ...] with header.")
                         }
                         init_cols <- seq(3,ncol(rawData),2)
                         dep_cols <- seq(4, ncol(rawData), 2)
                         init_counts <- rawData[, init_cols, with=F]
                         dep_counts <- rawData[, dep_cols, with=F]
                       } else {
                         init_counts <- NA
                         dep_counts <- rawData[, c(-1,-2), with=F]
                       }
                       if (any(duplicated(names(dep_counts)))) {
                         stop("please provide unique sample names.")
                       }

                       num_samples <- ncol(dep_counts)
                       
                       
                       # Read master data.
                       # use ModelObj to combine replicates etc.
                       private$write_log('Reading master data.')
                       if (any(is.na(masterFiles)) & !hasInitSeq) {
                         private$write_log("Must provide initial counts or master library abundances.")
                         stop("Must provide initial counts or master library abundances.")
                       }
                       if (any(is.na(masterFiles))) {
                         master_counts <- NA
                         if (length(masterFiles)==1) {
                           private$write_log('No master library provided.')
                         } else {
                           private$write_log("Invalid master library files included.")
                           private$write_log(masterFiles)
                           private$write_log('Warning: Invalid master lib files; not used.')
                           warning("Invalid master library files included; no master library used.")
                         }
                       } else {
                         private$write_log("submitted master data files are:")
                         private$write_log(masterFiles)
                         private$write_log('If duplicate guides found, both discarded.')
                         master_counts <- lapply(seq_along(masterFiles),
                                                 function(i){
                                                   rawdt <- fread(masterFiles[[i]])
                                                   dt <- unique(rawdt, by=1)
                                                   setnames(dt, old = c(1),
                                                            new = c('sgrna'))
                                                   return(dt)})
                         mf_names <- lapply(masterFiles,
                                            function(i) strsplit(i, '/')[[1]])
                         names(master_counts) <- sapply(mf_names,
                                                        function(i) strsplit(i[length(i)], '[.]')[[1]][1])
                         mf_illegal <- sapply(master_counts, function(i) {
                           !all(sapply(i[, -1, with=F], is.numeric))  |
                             any(sapply(i[, -1, with=F], function(j) j < 0)) 
                         })
                         if (any(mf_illegal)) {
                           private$write_log(lapply(master_counts[mf_illegal], head))
                           private$write_log("Invalid count data in masterlibrary counts!")
                           private$write_log('Error: masterFile must have form [sgrna, count1, count2..]')
                           stop("masterFile must have form [sgrna, count1, count2..]")
                         }
                       }
                       self$master_counts <- master_counts
                       
                       # make column lookup table for sample-masterlib mapping
                       if (any(is.na(masterFiles))) {
                         info_file <- NA
                       } else if (!is.na(sampleMasterInfoFile)) {
                         info_file <- fread(sampleMasterInfoFile, header = T)
                         if (!ncol(info_file) %in%  c(2,3)) {
                           stop("please submit sampleMasterInfoFile as 2 or 3 columns w/ header:
                                            sample, masterlib\nor\nsample, sampleType, masterlib")
                         }
                         if (ncol(info_file)==2) {
                           names(info_file) <- c("sample", "masterlib")
                         } else {
                           names(info_file) <- c('sample', 'sampleType', 'masterlib')
                         }

                         # if masterlib names entered with file extension, remove it.
                         mf_names <- lapply(info_file$masterlib,
                                            function(i) strsplit(i, '/')[[1]])
                         info_file[, masterlib := sapply(mf_names, function(i) {
                           strsplit(i[length(i)], '[.]')[[1]][1]})]
                         if (!all(info_file$masterlib %in% names(self$master_counts))) {
                           stop('Master count file names can not be matched to sampleMasterInfoFile.')
                         }
                         setkey(info_file, sample)
                       } else if (length(masterFiles) == 1) {
                         sampleNames <- names(dep_counts)
                         info_file <- data.table('sample' = sampleNames,
                                                 'masterlib' = names(master_counts))
                         private$write_log('Using one master library for all samples.')
                         private$write_log(head(info_file))
                         private$write_log(head(sampleNames))
                         setkey(info_file, sample)
                       } else {
                         stop("please provide assignment of master libraries to samples")
                       }
                       
                       # Use only count samples indicated in useSamples argument.
                       # meant to be an integer for single-sample analysis.
                       if (any(is.na(useSamples))) {
                         useSamples <- 1:ncol(dep_counts)
                       } else if (!(is.vector(useSamples, mode = 'integer') | is.integer(useSamples))) {
                         stop('supply useSamples as integer indices.')
                       } else {
                         # Select samples using indices relative to original count file. 
                         if (max(useSamples) > ncol(dupRawData)) {
                           private$write_log('Indices in useSamples do not correspond to count file.')
                           stop('Indices in useSamples do not correspond to count file.')
                         } else {
                           useSampleNames <- names(dupRawData)[useSamples]
                           useSampleCols <- which(names(dep_counts) %in% useSampleNames)
                         }
                       }
                       
                       
                       if (is.data.table(init_counts)) {
                         self$init_counts <- init_counts[, (useSampleCols), with=F]
                       } else {
                         self$init_counts <- NA
                       }

                       self$dep_counts <- dep_counts[, (useSampleCols), with=F]
                       
                       if (!is.data.table(info_file)) self$sample_masterlib <- NA
                       else {
                         self$sample_masterlib <- info_file
                       }
                       
                       # make sure all samples appropriately annotated, if relevant.
                       # take intersect of annotated, masterlib-annotated, and count-data-included
                       # sample names, and subset all data accordingly.
                       finalSampleNames <- names(self$dep_counts)
                       if (is.data.table(self$sample_annotations)) {
                         finalSampleNames <- intersect(finalSampleNames, self$sample_annotations$sample_name)
                       }
                       if (is.data.table(self$sample_masterlib)) {
                         finalSampleNames <- intersect(finalSampleNames, self$sample_masterlib$sample)
                       }
                       private$write_log(c('using samples: ', finalSampleNames))
                       message(c('Using ', length(finalSampleNames), ' samples.'))
                       
                       if (is.data.table(self$sample_annotations)) {
                         self$sample_annotations <- self$sample_annotations[sample_name %in% finalSampleNames, .SD]
                       } else { 
                         private$write_log('No annotation file used.')
                       }
                       if (is.data.table(self$sample_masterlib)) {
                         self$sample_masterlib <- self$sample_masterlib[sample %in% finalSampleNames, .SD]
                       } else {
                         private$write_log('No master library used.')
                       }
                       if (is.data.table(self$init_counts)) {
                         finalInitNames <- which(names(self$dep_counts) %in% finalSampleNames)
                         self$init_counts[, (finalInitNames), with=F]
                       } else {
                         private$write_log('No initial sequencing counts used.')
                       }
                       self$dep_counts <- self$dep_counts[, .SD, .SDcols = finalSampleNames]

                       
                       
                       if (is.na(guideCovarFile)) {
                         private$write_log("Not calculating guide efficiency. No covariate file.")
                         message("Not calculating guide efficiency. No covariate file.")
                       } else {
                         self$guide_covars <- guide_covars
                       }
                       
                       # Load guide2gene_map
                       # note guide and gene order MUST BE PRESERVED
                       # from rawData; used as lookup guide.
                       private$write_log("Using these columns as sgrna and gene labels:")
                       private$write_log(names(rawData)[1:2])
                       message(c("Using these columns as sgrna and gene labels:\n",
                               names(rawData)[1],'\t',names(rawData)[2]))
                       self$guide2gene_map <- rawData[,c(1,2)]
                       names(self$guide2gene_map) <- c("sgrna", "gene")
                       g_idx <- data.table("uniq_gene"=unique(self$guide2gene_map$gene))
                       g_idx[,gene_idx := 1:nrow(g_idx)]
                       setkey(g_idx, uniq_gene)
                       self$guide2gene_map[, gene_idx := g_idx[gene, gene_idx]]
                       
                       user_total_reads = as.numeric(totalReads)
                       if (!any(is.na(user_total_reads))) {
                         num_samples_given <- length(user_total_reads)
                         if (num_samples_given == 1) {
                           self$init_total_reads <- rep(user_total_reads, num_samples)
                           self$dep_total_reads <- rep(user_total_reads, num_samples)
                         } else if (num_samples_given == 2) {
                           self$init_total_reads <- rep(user_total_reads[[1]], num_samples)
                           self$dep_total_reads <- rep(user_total_reads[[2]], num_samples)
                         } else if (num_samples_given == num_samples) {
                           self$init_total_reads <- user_total_reads
                           self$dep_total_reads <- user_total_reads
                         } else if (num_samples_given == 2 * num_samples) {
                           private$write_log("Assuming total reads submitted as pairs of initial and depleted.")
                           message("Assuming total reads arg submitted as pairs of initial and depleted.")
                           self$init_total_reads <- user_total_reads[seq(1, num_samples_given, 2)]
                           self$dep_total_reads <- user_total_reads[seq(2, num_samples_given,2)]
                         } else {
                           private$write_log(c('Error: the number of total_reads ',
                                               "does not make sense for this many samples:",
                                               num_samples_given))
                           stop(c("The number of total_reads does not make sense for this many samples:",
                                  num_samples_given))
                         }
                       } else if (hasInitSeq) {
                         self$init_total_reads <- colSums(self$init_counts)
                         self$dep_total_reads <- colSums(self$dep_counts)
                       } else {
                         self$init_total_reads <- NA
                         self$dep_total_reads <- colSums(self$dep_counts)
                       }
                       
                       if (is.na(cellsInfected)) {
                         self$cells_infected <- rep(1e3*nrow(rawData), num_samples)
                       } else if (length(cellsInfected==1)) {
                         self$cells_infected <- rep(cellsInfected * nrow(rawData), num_samples)
                       } else {
                         self$cells_infected <- cellsInfected
                       }
                       
                       # check for 0-counts across samples, guides, & genes,
                       # and remove data accordingly.
                       # Add pseudocount, & reset keys afterward.
                       . <- removeNullData(self, private$write_log)
                       
                       # Read in file of negative control genes.
                       # Used in calculating sample scaling factor
                       if (!is.na(negCtrlFile)) {
                         private$write_log(c('reading negative control file.',
                                             '\nFile sample (should have header):',
                                             head(negCtrlFile)))
                         self$neg_control_file <- negCtrlFile
                       } else self$neg_control_file <- NA
                       message('DataObj initialized.')
                     }
                   )
)
# -------------------------------- End of Script ---------------------------------

