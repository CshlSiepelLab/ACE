#' DataObjClass: R6 Environment Object for Raw Data
#'
#' This function is an R6 object class to create an environment for data analysis.
#' Use to load all data associated with a single-target negative selection
#' CRISPR screen; stores all data with minimal processing, but checks for
#' errors and data inconsistencies.
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


                                 #' @description Remind users of arguments.
                                 help = function() {
                                   print("#' newObj <- DataObj$new(masterFile.txt, countFile.txt,
                                cells_infected, [guideCovarFile.txt],
                                [sampleInfo.txt], [totalReads]) ")
                                 },

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
                                 #' @param useSamples Sample indices to analyze.
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

                                   # Debugging.
                                   if (R.version$major < 3 || (R.version$major == 3 & R.version$minor < 3)) {
                                     warning("R version too old; please load version 3.3.0 or higher.")
                                   }
                                   if (packageVersion('data.table')$minor < 11) {
                                     print("package version of data.table is:")
                                     print(packageVersion('data.table'))
                                     warning("Please install the latest version of data.table.")
                                   }

                                   # Read optional guide and sample information files
                                   if (!is.na(guideCovarFile)) {
                                     guide_covars <- fread(guideCovarFile)
                                     if (ncol(guide_covars) != 2) {
                                       print("submitted sgrna covariate file is:")
                                       print(head(guide_covars))
                                       stop("seqFile must have form [sgrna, sequence] with header;
                                            may add more covariates later.")
                                     }
                                   }
                                   # have semicolon-delimited 2nd column, with sample anntations.
                                   if (!is.na(sampleInfoFile)) {
                                     rawSampleInfo <- tryCatch(fread(sampleInfoFile),
                                                               warning = function(i) simpleError('Invalid sampleInfoFile'))

                                     if (ncol(rawSampleInfo) != 2) {
                                       print("submitted sample metadata file is:")
                                       print(head(rawSampleInfo))
                                       stop("sample meta data must be in form [sample_name, sample_annotations]")
                                     }
                                     names(rawSampleInfo) <- c("sample_name", "sample_subtype")
                                     self$sample_annotations <- rawSampleInfo
                                   }

                                   # Load Count Data
                                   dupRawData <- tryCatch(fread(countFile, header = T),
                                                          warning = function(i) simpleError('Invalid countData file'))
                                   # check for duplicates, keep first instance.
                                   if (any(duplicated(dupRawData[[1]]))) {
                                     print('using first column as UNIQUE guide labels.')
                                     print(head(dupRawData))
                                     print('Some guides listed twice; discarding')
                                     print(c('See dup_guides.txt for all duplicated counts;',
                                             'total are:'))
                                     print(sum(duplicated(dupRawData[[1]])))
                                     write.table(dupRawData[duplicated(dupRawData[[1]]),],
                                                 file = 'dup_guides.txt',
                                                 quote = F)
                                     rawData <- unique(dupRawData, by = 1)
                                   } else {
                                     rawData <- dupRawData
                                   }
                                   setorderv(rawData, names(rawData)[2])

                                   # check for non-numeric data.
                                   if (any(sapply(rawData[,-(1:2), with=F],
                                                   function(i) !is.na(i) & !is.numeric(i)))) {
                                     write.table(rawData[sapply(.SD, !is.na) &
                                                           sapply(.SD, !is.numeric),
                                                         .SDcols = -(1:2)],
                                                 quote = F,
                                                 row.names = F,
                                                 file = 'data_debug.txt')
                                     stop("non-numeric and not-NA data in raw count file.")
                                   }
                                   # split into initial & dep counts if applicable.
                                   if (hasInitSeq) {
                                     if (ncol(rawData) < 4 || ncol(rawData)%%2 == 1) {
                                       print("input count file is:")
                                       print(head(rawData))
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
                                   if (any(is.na(masterFiles)) & !hasInitSeq) {
                                     stop("Must provide initial counts or master library abundances.")
                                   }
                                   if (any(is.na(masterFiles))) {
                                     master_counts <- NA
                                     print("No list of master library files included.")
                                   } else {
                                     print("submitted master data files are:")
                                     print(masterFiles)
                                     print('If duplicate guides found, first instance used.')
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
                                       !all(sapply(i[, -1, with=F], is.numeric)) })
                                     if (any(mf_illegal)) {
                                       print(lapply(master_counts[mf_illegal], head))
                                       print("Non-numeric in masterlibrary counts!")
                                       stop("masterFile must have form [sgrna, count1, count2..]")
                                     }
                                   }
                                   self$master_counts <- master_counts

                                   # make column lookup table for samples
                                   sampleNames <- names(dep_counts)
                                   if (any(is.na(masterFiles))) {
                                     info_file <- NA
                                     print(head(info_file))
                                     # TODO: BUG from line 256 leads here.
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
                                     if (!all(sampleNames %in% info_file$sample)) {
                                       print(sampleNames[which(!sampleNames %in% info_file$sample)])
                                       print(head(info_file))
                                       stop("sampleMasterInfoFile needs entry per count sample")
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
                                   } else if (length(masterFiles == 1)) {
                                     info_file <- data.table('sample' = sampleNames,
                                                             'masterlib' = names(master_counts))
                                     print(head(info_file))
                                     print(sampleNames)
                                     setkey(info_file, sample)
                                   } else {
                                     stop("please provide assignment of master libraries to samples")
                                   }

                                   # Threshold data according to useSamples argument.
                                   # meant to be an integer for single-sample analysis.
                                   # TODO: organize all masterlib~countdata sample name checking into one section.
                                   if (is.na(useSamples)) useSamples <- 1:ncol(dep_counts)
                                   if (!(is.vector(useSamples, mode = 'integer') | is.integer(useSamples))) {
                                     stop('supply useSamples as integer indices.')
                                   }
                                   if ('sampleType' %in% info_file) {
                                     useSampleNames <- info_file[sampleType %in% unique(sampleType)[useSamples], sample]
                                     useSampleCols <- which(names(dep_counts) %in% useSampleNames)
                                   } else {
                                     useSampleCols <- useSamples
                                   }
                                   if (is.data.table(init_counts)) {
                                     self$init_counts <- init_counts[, (useSampleCols), with=F]
                                   } else {
                                     self$init_counts <- NA
                                   }
                                   self$dep_counts <- dep_counts[, (useSampleCols), with=F]
                                   cat('using samples: ', names(dep_counts[useSampleCols]))
                                   finalSampleNames <- names(self$dep_counts)
                                   
                                   if (!is.data.table(info_file)) self$sample_masterlib <- NA
                                   else {
                                     self$sample_masterlib <- info_file[sample %in% finalSampleNames,]
                                   }

                                         
                                   

                                   if (is.na(guideCovarFile)) {
                                     print("Not calculating guide efficiency. No covariate file.")
                                   } else {
                                     self$guide_covars <- guide_covars
                                   }

                                   # Load guide2gene_map
                                   # note guide and gene order MUST BE PRESERVED
                                   # from rawData; used as lookup guide.
                                   print("Using these columns as sgrna and gene labels.")
                                   print(names(rawData)[1:2])
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
                                       print("Assuming total reads submitted as pairs of initial and depleted.")
                                       self$init_total_reads <- user_total_reads[seq(1, num_samples_given, 2)]
                                       self$dep_total_reads <- user_total_reads[seq(2, num_samples_given,2)]
                                     } else {
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
                                   . <- removeNullData(self)

                                   # Read in file of negative control genes.
                                   # Used in calculating sample scaling factor
                                   if (!is.na(negCtrlFile)) {
                                     print('reading negative control file; please include header.')
                                     self$neg_control_file <- negCtrlFile
                                   } else self$neg_control_file <- NA
                                 }
                   )
)
# -------------------------------- End of Script ---------------------------------

