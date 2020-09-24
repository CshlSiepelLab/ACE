#' Function ResultObjClass
#'
#' An R6 object class to present the results of model fitting as produced by optimizeModelParameters.
#' Not called by user, but rather within optimizeModelParameters.
#' @docType class
#' @importFrom R6 R6Class
#' @import data.table
#' @param optObjList List from optim fit. With null_liklihood, etc. added on.
#' @param DataObj Produced from raw data.  just need guide2gene_map and num_samples though.
#' @param ModelObj Inferred parameteer object, just need mean_var_model.
#' @return ResultObj R6 class object, with gene_essentiality, guide_efficiency, sample_response,
#' guide_covariates, gene_essentiality_likelihood, gene_essentiality_CI, gene_essentiality_pvalue,
#' differential_depletion_pvals, mean_var_model_parameters.
#' @export
ResultObj <- R6Class("ResultObj",
                     public = list(
                       #' @field gene_essentiality Data.table of gene essentiality
                       #' for all, test, and ctrl sample panels.
                       gene_essentiality = "data.table",
                       #' @field guide_efficiency Vector of sgRNA feature names used.
                                   guide_efficiency = "vector",
                       #' @field guide_covariate_score Matrix of sgRNA covariates.
                                   guide_covariate_score = "matrix",
                       #' @field guide_covariate_weight Vector of sgRNA feature weights.
                                   guide_covariate_weight = "vector",
                       #' @field sample_response Vector of inferred sample scaling factors.
                                   sample_response = "vector",
                       #' @field gene_essentiality_likelihood Data.table of likelihood values
                       #' for all, all_null, test, test_null, ctrl, and ctrl_null tests.
                                   gene_essentiality_likelihood = "data.table",
                       #' @field gene_essentiality_low_CI Lower confidence interval of gene essentiality estimates in \code{gene_essentiality}.
                                   gene_essentiality_low_CI = "data.table",
                       #' @field gene_essentiality_hi_CI Higher confidence interval of gene essentiality estimates in \code{gene_essentiality}.
                                   gene_essentiality_hi_CI = "data.table",
                       #' @field gene_essentiality_pvalue Empirical p-values for within- and between-panel essentiality; data.table.
                                   gene_essentiality_pvalue = "data.table",
                       #' @field mean_var_model_parameters Vector of inferred mean-variance relationship parameters; model specific.
                                   mean_var_model_parameters = "vector",

                       #' @description Create ResultObj.
                     initialize = function(optObjList, DataObj, ModelObj) {
                       # include essentiality == 1 in pval list? Not required for diff depletion.
                       # 0 or 1 value.
                       gene_names <- DataObj$guide2gene_map
                       num_genes <- length(gene_names)
                       num_samples <- ncol(DataObj$init_counts)
                       num_guide_covars <- length(optObjList[["all"]]$par) - num_genes - num_samples

                       # load as genes x [all, test, ctrl] data.table
                       self$gene_essentiality <- data.table("Gene" = gene_names,
                                                            "All_Samples" = optObjList[["all"]]$par[1:num_genes],
                                                            "Test_Samples"= optObjList[["test"]]$par[1:num_genes],
                                                            "Ctrl_Samples"= optObjList[["ctrl"]]$par[1:num_genes])
                       self$guide_efficiency <- optObjList[["all"]]$guide_efficiency
                       self$guide_covariate_score <- optObjList[["all"]]$guide_covariate_score
                       self$guide_covariate_weight <- optObjList[["all"]]$par[(num_genes+1):(num_genes+1+num_guide_covars)]

                       self$sample_response <- optObjList[["all"]]$par[-c(1:(num_genes+num_guide_covars+2))]

                       self$gene_essentiality_likelihood <- data.table("Gene" = gene_names,
                                                                        "All_Samples" = optObjList[["all"]]$value,
                                                                        "Test_Samples"= optObjList[["test"]]$value,
                                                                        "Ctrl_Samples"= optObjList[["ctrl"]]$value)
                       self$guide_null_likelihood <- data.table("Gene" = gene_names,
                                                                "All_Samples" = optObjList[["all"]]$null_likelihood,
                                                                "Test_Samples"= optObjList[["test"]]$null_likelihood,
                                                                "Ctrl_Samples"= optObjList[["ctrl"]]$null_likelihood)
                       getCI <- function(hessian_obj) {
                         ste <- -sqrt(-hessian_obj)
                         return(1.96*ste)
                       }
                       getPval <- function(logAlt, logNull, df) {
                         pval <- 1 - pchisq(2 * (logAlt - logNull), df)
                       }

                       #TODO: make sure hessian can be subset for JUST the gene parameters.
                       all_CI_95 <- getCI(optObjList[["all"]]$hessian[1:num_genes])
                       test_CI_95<- getCI(optObjList[["test"]]$hessian[1:num_genes])
                       ctrl_CI_95 <- getCI(optObjList[["ctrl"]]$hessian[1:num_genes])
                       all_low <- self$gene_essentiality$All_Samples - all_CI_95
                       all_hi <- self$gene_essentiality$All_Samples + all_CI_95
                       test_low <- self$gene_essentiality$Test_Samples - test_CI_95
                       test_hi <- self$gene_essentiality$Test_Samples + test_CI_95
                       ctrl_low <- self$gene_essentiality$Ctrl_Samples - test_CI_95
                       ctrl_hi <- self$gene_essentiality$Ctrl_Samples + test_CI_95
                       self$gene_essentiality_low_CI <- data.table("Gene" = gene_names,
                                                                   "All_Samples" = all_low,
                                                                   "Test_Samples" = test_low,
                                                                   "Ctrl_Samples" = ctrl_low)
                       self$gene_essentiality_hi_CI <- data.table("Gene" = gene_names,
                                                                  "All_Samples" = all_hi,
                                                                  "Test_Samples" = test_hi,
                                                                  "Ctrl_Samples" = ctrl_hi)
                        logAlt <- self$guide_essentiality_likelihood[,Test_Samples * Ctrl_Samples]
                       logNull <- self$guide_essentiality_likelihoood[, All_Samples]
                       self$differential_depletion_pvals <- getPval(logAlt, logNull, df=1)

                       all_pval <- getPval(logAlt = self$gene_essentiality_likelihood$All_Samples,
                                           optObjList[["all"]]$null_likelihood, df = 1)
                       test_pval <- getPval(logAlt = self$gene_essentiality_likelihood$Test_Samples,
                                            optObjList[["test"]]$null_likelihood, df = 1)
                       ctrl_pval <- getPval(logAlt = self$gene_essentiality_likelihood$Ctrl_Samples,
                                            optObjList[["ctrl"]]$null_likelihood, df = 1)
                       self$gene_essentiality_pvalue <- data.table("Gene" = gene_names,
                                                                   "All_Samples"=all_pval,
                                                                   "Test_Samples"=test_pval,
                                                                   "Ctrl_Samples"=ctrl_pval)
                       self$mean_var_model_parameters <- ModelObj$mean_var_model
                     }
                     )
                     )
