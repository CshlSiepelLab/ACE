#' ModelObj Environment
#'
#' @description 
#' This R6 environment object processes the raw data loaded in the provided
#' \code{'user_DataObj'} to calculate model constants according to user
#' specifications.  Computes sample scaling parameters, the vector of infected
#'  cell values to use in the likelihood summation,
#'  and partitions test and control samples,
#' @docType class
#' @import R6
#' @importFrom R6 R6Class
#' @import data.table
#' @importFrom stringr str_count
#' @include callGetLLByGene.R
#' @param user_DataObj Containing all count data; a DataObj.
#' @param use_master_library Boolean, use masterlibrary to set guide abundance prior?
#' @param fit_guide_parameter Boolean, should by-guide efficiencies be calculated?
#' @param mean_var_model Optional integer for mean~var fitting model; default poisson.
#' @param use_neg_ctrl Boolean, use negative control genes for scaling sample
#'          parameters and final essentiality value.
#' @param test_samples String indicating which sample annotation to use when
#'          grouping test and control groups.
#' @param stepSize Number of cells to use for each step in the summation, default 10.
#' @return ModelObj (R6Class).
#' Contains:
#' --  mean_var_model_params = "list",
#' --  guide_features = "matrix",
#' --  unobserved_infected_cell_values = "vector",
#' --  mean_var_model = "integer",
#' --  init_scaling = "vector",
#' --  dep_scaling = "vector",
#' --  test_sample_subtype_cols = "vector"
#' -- stepSize = 10
#' -- guidePrior = 'string'
#' @export
ModelObj <- R6Class("ModelObj",
                    private = list(
                      write_log = 'function'),
                    public = list(
                      #' @field mean_var_model_params Coefficients to use in mean-var model.
                      mean_var_model_params = "list",
                      #' @field guide_features Matrix of by-guide features.
                      guide_features = "matrix",
                      #' @field master_freq_dt Data.table of sgRNA frequencies in master libraries.
                      master_freq_dt = 'data.table',
                      #' @field unobserved_infected_cell_values Discrete values of infected cells for Riemann sum.
                      unobserved_infected_cell_values = "vector",
                      #' @field mean_var_model Integer specifying which form of model to use.
                      mean_var_model = "integer",
                      #' @field init_scaling By-sample vector of sample scaling factors.
                      init_scaling = "vector",
                      #' @field dep_scaling By-sample vector of depleted sample scaling.
                      dep_scaling = "vector",
                      #' @field test_sample_subtype_cols Vector of column indices for 'test' subtype.
                      test_sample_subtype_cols = 'vector',
                      #' @field stepSize Integer of infected cell intervals to evaluate.
                      stepSize = 10,
                      #' @field guidePrior String indicating form of infection prior for README.
                      guidePrior = 'string',
                      #' @field use_neg_ctrl Boolean indicating whether negative controls from DataObj should be used.
                      use_neg_ctrl = T,
                      #' @field neg_ctrls String vector of gene names for negative controls.
                      neg_ctrls = 'vector',

                      #' @description Create ModelObj
                      initialize = function(user_DataObj=NA,
                                            use_master_library = T,
                                            fit_guide_parameter = F,
                                            mean_var_model = 1,
                                            use_neg_ctrl=T,
                                            test_samples = NA,
                                            stepSize=10) {
                        
                        # Create log file for info, warnings, and errors.
                        if (!dir.exists('ACE_output_data')) dir.create('ACE_output_data')
                        tStamp <- paste(unlist(str_split(Sys.time(), ' |:')), collapse='_')
                        log_file <- file(file.path('ACE_output_data',
                                                   paste0('ACE_ModelObj_log_', tStamp,
                                                                            '.txt.')),
                                         open='w+')
                        on.exit(close(log_file),add=T)
                        
                        #'  Write messages and warnings to log file.
                        #'  @param message_vector String or table to write.
                        private$write_log <- function(message_vector) {
                          if (is.atomic(message_vector)) {
                            cat(message_vector,'\n', file=log_file, append=T)
                          } else {
                            suppressWarnings(write.table(message_vector, file = log_file, 
                                        col.names=T, append=T))
                          }
                        }
                        private$write_log('Log File for ModelObj Generation')
                        
                        if (!is.environment(user_DataObj)) {
                          message('initialize DataObj with:\n',
                                  "ModelObj$new(user_DataObj, use_master_library=T, ",
                                  "fit_guide_parameter=F, mean_var_model=1, use_neg_ctrl=T, ",
                                  "test_samples = NA)")
                          return()
                        }

                        # fit mean-var model
                        # normalize init_counts by method of means (median of ratio of counts/mean)
                        user_mean_var_model <- mean_var_model
                        if (!is.data.table(user_DataObj$init_counts)) {
                          private$write_log('Using masterlibrary and depletion counts.')
                          use_counts <- user_DataObj$dep_counts
                        } else {
                          use_counts <- user_DataObj$init_counts
                        }
                        if (user_mean_var_model ==0) {
                          self$mean_var_model <- 0L
                          self$mean_var_model_params <- getMeanVarFit0(use_counts)
                        } else if (user_mean_var_model == 1 || length(grep(user_mean_var_model, "Poisson",
                                                                           ignore.case=T))) {
                          private$write_log("Using poisson distribution for sequencing noise.")
                          self$mean_var_model <- 1L
                          self$mean_var_model_params <- list(0, 0)
                        } else if (user_mean_var_model == 2 || length(grep(user_mean_var_model,
                                                                           "trend",
                                                                           ignore.case=T))) {
                          private$write_log("Using glm fit mean~var model for NB distribution.")
                          self$mean_var_model <- 2L
                          self$mean_var_model_params <- getMeanVarFit2(use_counts)
                        } else if (user_mean_var_model == 3 || length(grep(user_mean_var_model,
                                                                           c("ribodiff", "bayes","NB"),
                                                                           ignore.case))) {
                          private$write_log("Using arg max mean~var fit model for NB distribution.")
                          self$mean_var_model <- 3L
                          self$mean_var_model_params <- getMeanVarFit3(use_counts)
                        } else if (user_mean_var_model == 4 || length(grep(user_mean_var_model, "DESEQ",
                                                                           ignore.case=T))) {
                          private$write_log("Using DESEQ's variance model.")
                          self$mean_var_model <- 4L
                          self$mean_var_model_params <- getMeanVarFit4(use_counts)
                        } else {
                          private$write_log('ERROR: invalid mean~var model specified.')
                          stop("Invalid mean~var model specified.")
                        }

                        # Normalize guide features to same mean and range (-1/2 to 1/2).
                        # Feature weights with sigmoid transform will be used to
                        # convert to percent of 'efficient' guides.
                        # Can be equivalently reasoned as the fraction of 
                        # effectively cutting guides or the penetrance of the
                        # gene fitness effect.
                        # first column of user_DAtaObj$guide_covar is sgRNA ID.
                        if (fit_guide_parameter & is.data.table(user_DataObj$guide_covar)) {
                          guide_features <- apply(user_DataObj$guide_covar[,.SD, 
                                                                           .SDcols=-1],
                                                  2, function(i) {
                                                    (i - mean(i))/(max(i)-min(i))
                                                  })
                          isVaried <- apply(guide_features, 2, uniqueN)
                          if (all(isVaried < 2)) {
                            message("all guides have an identical covariate; unable to fit guide parameter.")
                            private$write_log("all guides have identical covariates; unable to fit guide parameter.")
                            guide_features <- rep(1, nrow(user_DataObj$guide_covar))
                          } else if (any(isVaried < 2)) {
                            message('Some covariates identical across all guides; discarding.')
                            private$write_log('Some covariates identical across all guides; discarding.')
                            guide_features[, .SD := NULL, .SDcols = isVaried < 2]
                          }
                          guide_features[, 'sgrnaID':= userDataObj$guide_covar[[1]]]
                          setcolorder(guide_features, 'sgrnaID')
                          self$guide_features <- guide_features[match(sgrnaID,
                                                                      user_DataObj$guide2gene_map$sgrna),]
                          # TODO: check for categorical features & expand.
                          # TODO: check for NA's in submitted file.
                          message("Using the following features to fit guide efficiency:",
                                  names(self$guide_features)[-1])
                          private$write_log("Using the following features to fit guide efficiency:",
                                  names(self$guide_features))
                          
                        } else {
                          self$guide_features <- NA
                        }

                        # n_sg values; represent number of infected cells.
                        # Assume max single guide enrichment is min(10%, 100*avg_cells)
                        # and use a resolution of 10 cells.
                        self$stepSize <- stepSize
                        maxCells <- min(100*user_DataObj$cells_infected[[1]]/nrow(use_counts),
                                        round(0.1 * user_DataObj$cells_infected[[1]]))
                        self$unobserved_infected_cell_values <- seq(stepSize, maxCells, stepSize)
                        message('Summation performed over vector of length ',
                                length(self$unobserved_infected_cell_values), '\n')
                        private$write_log(c('Summation performed over vector of length ',
                                length(self$unobserved_infected_cell_values), '\n'))

                        # get vector of 'test' columns for convenience.
                        # using keywords in second column of sample annotation file.
                        if (is.na(test_samples) || !is.data.table(user_DataObj$sample_annotations)) {
                          message("Not calculating differential depletion by sample subtype.")
                          private$write_log("Not calculating differential depletion by sample subtype.")
                        } else {
                          if (length(test_samples) > 1) {
                            message('One testing group per model;',
                                                            'using first listed feature only.\n')
                            private$write_log('One testing group per model;',
                                                            'using first listed feature only.\n')
                          }
                          message(c("extracting sample subtype:", test_samples[1]))
                          private$write_log(c("extracting sample subtype:", test_samples[1]))
                          
                          test_sample_rows <- grep(test_samples[1],
                                                   user_DataObj$sample_annotations$sample_subtype)
                          test_sample_names <- user_DataObj$sample_annotations$sample_name[test_sample_rows]
                          private$write_log("test samples in sample annotations file:")
                          private$write_log(test_sample_names)

                          # get list of test samples; must correspond
                          # to sample_masterlib listings, which are
                          # depleted sample names.
                          getSampleIdx <- function(sampleName) {
                            ifelse(sampleName %in% names(user_DataObj$dep_counts),
                                   grep(sampleName, names(user_DataObj$dep_counts)),
                                   grep(sampleName, names(user_DataObj$init_counts)))
                          }
                          anno_sample_idx <- sapply(test_sample_names, getSampleIdx)
                          test_sample_idx <- unique(anno_sample_idx[!is.na(anno_sample_idx)])

                          if (length(test_sample_idx)==0) {
                            message("No samples matching the test subtype found in count data;",
                                              "no differential essentiality will be calculated.")
                            private$write_log(c("No samples matching the test subtype found in count data;",
                                              "no differential essentiality calculated."))
                          } else if (length(test_sample_idx)==ncol(user_DataObj$dep_counts)) {
                            message(paste0('All samples in count data match the test subtype;',
                                           ' no differential essentiality calculated.'))
                            private$write_log(paste0('All samples in count data ',
                                                     'match the test subtype;',  
                                                     'no differential essentiality calculated.'))
                          } else {
                            self$test_sample_subtype_cols <- test_sample_idx
                            message("Using the following number of samples as our test set:")
                            message(length(self$test_sample_subtype_cols))
                            private$write_log("Using the following samples as our test set:")
                            private$write_log(c(names(user_DataObj$dep_counts)[self$test_sample_subtype_cols]))
                          }
                        }

                        # n_sg prior reflects:
                        # - masterlib * # cells infected
                        # - uniform across all values of n_sg (1/c_s)
                        # - by-guide median abundance * # cells infected.
                        # n_sg will be scaled to median of read counts using γ
                        #  If master counts present:
                        # Combine replicates of each masterlib into a single by-guide
                        # vector of frequencies.
                        # when calling multiple samples, entire table must be passed.
                        if (is.data.table(user_DataObj$master_counts[[1]]) & use_master_library) {
                          private$write_log("Using master library for 'cells infected' prior calculation")
                          # get data.table of sgrna frequencies averaged across masterlibs.
                          master_freq_dt <- getMasterLib(countList = user_DataObj$master_counts,
                                                         private$write_log)
                          self$master_freq_dt <- master_freq_dt
                          self$guidePrior <- 'master_library'
                          private$write_log('Model Obj master_freq_dt is:')
                          private$write_log(head(self$master_freq_dt))
                        } else if (is.data.table(user_DataObj$init_counts)) {
                          private$write_log("Using median abundance in initial counts as prior for cells infected.")
                          # assume no counts are depleted in the initial population.
                          guideMeans <- user_DataObj$init_counts[, exp(rowMeans(log(.SD+0.5)))]
                          normCounts <- user_DataObj$init_counts[,(0.5+.SD)/
                                                                   sapply(seq_along(.SD),
                                                                          function(i)
                                                                            median((0.5+user_DataObj$init_counts[[i]])/
                                                                                     guideMeans[i]))]
                          self$master_freq_dt<-data.table('sgrna' = user_DataObj$guide2gene_map$sgrna,
                                                          'meanFreq' = rowMeans(
                                                            normCounts[,log(.SD)]-
                                                              user_DataObj$init_counts[,log(colSums(0.5+.SD))]))
                          setkey(self$master_freq_dt, sgrna)
                          self$guidePrior <- 'mean_initial'
                        } else {
                          private$write_log('Using uniform prior for cells infected.')
                          self$master_freq_dt <- NA
                          self$guidePrior <- 'uniform_cells'
                        }

                        # if negative controls present, use them.
                        if (is.character(user_DataObj$neg_control_file) & use_neg_ctrl) {
                          self$use_neg_ctrl <- T
                          # Remove outlying negative controls according to LFC (ACE agnostic).
                          # keep only controls with a LFC within 1 standard deviation of the median.
                          neg_ctrls <- fread(user_DataObj$neg_control_file, header =F,
                                             sep = '')[[1]]
                          private$write_log('treating negative control file as single column.')
                          private$write_log(head(neg_ctrls))
                          depSampleNames <- names(user_DataObj$dep_counts)
                          if (is.data.table(user_DataObj$init_counts)) {
                            use_base_counts <- user_DataObj$init_counts 
                          } else {
                            if (ncol(self$master_freq_dt)==2) {
                              useMasterSamples <- names(self$master_freq_dt)[2]
                            } else {
                              useMasterSamples <- sapply(depSampleNames, function(i) {
                                user_DataObj$sample_masterlib[sample == i, masterlib]
                              })
                              if (useMasterSamples[1] != 
                                user_DataObj$sample_masterlib[sample==depSampleNames[1], 
                                                              masterlib]){
                                stop('Incorrectly subset masterlibrary.')
                              }
                            }
                            use_base_counts <- self$master_freq_dt[
                              match(user_DataObj$guide2gene_map$sgrna, sgrna),.SD, 
                              .SDcols = useMasterSamples]                           
                          }
                          # get lfc for all genes.
                          avgNeg <- getLfcDt(user_DataObj = user_DataObj,
                                             getCI=F,
                                             isSim = F,
                                             use_base_counts = use_base_counts,
                                             write_log = private$write_log)
                          message('Removed negative controls more than one SD from mean LFC.')
                          sdNeg <- sd(avgNeg[gene %in% neg_ctrls, score])
                          medNeg <- median(avgNeg[gene %in% neg_ctrls, score])
                          useNeg <- avgNeg[gene %in% neg_ctrls &
                                             abs(score - medNeg) < sdNeg, gene]
                          self$neg_ctrls <- unique(useNeg)

                        } else {
                          message('No negative controls will be used.')
                          private$write_log('No negative controls will be used.')
                          self$use_neg_ctrl <- F
                          self$neg_ctrls <- NA
                        }

                        # Sample Scaling
                        sample_scaling <- deriveSampleScaling(user_DataObj,
                                                              self$master_freq_dt,
                                                              use_master_library,
                                                              self$use_neg_ctrl,
                                                              neg_ctrls = self$neg_ctrls,
                                                              private$write_log)
                        self$init_scaling <- sample_scaling$init_scaling
                        self$dep_scaling <- sample_scaling$dep_scaling

                        # for time saving measure, assume infected cells are within 3 orders of
                        # magnitude of the average.
                        if (max(user_DataObj$cells_infected) < nrow(use_counts)/100) {
                          print("cells used in the infection:")
                          print(max(user_DataObj$cells_infected))
                          print("guides infected:")
                          print(nrow(use_counts))
                          private$write_log(c('ERROR: "Not enough cells per guide;",
                                   "must infect at least 1/100 cells per guide in at least one sample\n'))
                          stop(cat("Not enough cells per guide;",
                                   "must have at least 1/100 cells per guide in at least one sample\n"))
                        }

                        # Debug everything calculated properly.
                        # mean_var_model_params can be a list of 2 numerics. Or a vector length guides.
                        # guide_features should be a numeric matrix [numGuides x numFeatures].
                        if (any(sapply(c(self$mean_var_model_params,
                                         self$unobserved_infected_cell_values, self$mean_var_model),
                                       function(i) any(is.na(unlist(i)))))) {
                          private$write_log("ERROR: NA's generated in model object in:")
                          if (any(is.na(self$mean_var_model))) private$write_log("mean_var_model")
                          if (any(is.na(self$mean_var_model_params))) private$write_log("mean_var_Model_params")
                          stop("Error in ModelObjClass; inappropriate NA's found in data.")
                        }
                      }
                      
                    )
)
# -------------------------------- End of Script -----------------------------------

