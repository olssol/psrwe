#' @title Outcome Analysis for PS-Integrated Estimation
#'
#' @description
#' Report outcome analysis for the PS-integrated approach.
#'
#' @param dta_psrst A returned object with class \code{PSRWE_EST}
#' @param method_ci A method name for confidence interval (default wald)
#' @param conf_type A type name of transformation for the confidence interval
#'        of PSKM approach (default log_log)
#' @param conf_int A two-sided level of confidence/credible limits
#'        (default 0.95)
#' @param alternative A character string for the alternative hypothesis that
#'        must be one of \code{"less"} (default) or \code{"greater"}, or
#'        \code{"two_sided"} (for log-rank and RMST only)
#' @param mu A number indicating the true value of the parameter of interest
#'        (or the difference in means for two arms),
#'        \code{mu = 0} when the test is log-rank or RMST
#' @param method_pval A method name for p-value (default wald),
#'        no impact for Bayesian method, and
#'        \code{method = "score"} only is for binary outcome in
#'        single arm study (i.e., comparing with a PG set by \code{mu})
#' @param ... Other options
#'
#' @return A list with class name \code{PSRWE_EST_OUTANA}.
#'
#' @details This function is mainly for summarizing and reporting the
#'     outcome analysis for the PS-integrated estimation.
#'     The input \code{dta_psrst} can be generated from the functions
#'     \code{\link{psrwe_powerp}}, \code{\link{psrwe_compl}}, and
#'     \code{\link{psrwe_survkm}}.
#'     See the functions \code{\link{psrwe_ci}} and \code{\link{psrwe_infer}}
#'     for the options of outcome analyses.
#'
#' @examples
#' data(ex_dta)
#' dta_ps <- psrwe_est(ex_dta,
#'        v_covs = paste("V", 1:7, sep = ""),
#'        v_grp = "Group",
#'        cur_grp_level = "current")
#' ps_borrow <- psrwe_borrow(total_borrow = 30, dta_ps)
#' ps_rst <- psrwe_compl(ps_borrow, v_outcome = "Y_Con")
#' rst <- psrwe_outana(ps_rst)
#' rst
#'
#' @export
#'
psrwe_outana <- function(dta_psrst,
                         method_ci = c("wald", "wilson"),
                         conf_type = c("log_log", "plain"),
                         conf_int = 0.95,
                         alternative = c("less", "greater", "two_sided"),
                         mu = 0,
                         method_pval = c("wald", "score", "score_weighted"),
                         ...) {
    ## check
    stopifnot(inherits(dta_psrst,
                       what = get_rwe_class("ANARST")))

    stopifnot(dta_psrst$Method %in% get_rwe_class("ANAMETHOD"))

    method_ci <- match.arg(method_ci)
    conf_type <- match.arg(conf_type)
    alternative <- match.arg(alternative)
    method_pval <- match.arg(method_pval)

    ## check components
    outcome_type <- dta_psrst$Outcome_type

    is_rct <- dta_psrst$is_rct
    if (is_rct) {
        type <- "Effect"
    } else {
        type <- "Control"
    }

    is_km <- dta_psrst$Method %in% get_rwe_class("ANAMETHOD_KM")
    col_est <- c("Mean", "StdErr")
    if (is_km) {
        col_est <- c(col_est, "T")
    }

    ## add ci and infer
    dta_psrst <- psrwe_ci(dta_psrst,
                          method_ci = method_ci,
                          conf_type = conf_type,
                          conf_int = conf_int,
                          ...)
    dta_psrst <- psrwe_infer(dta_psrst,
                             alternative = alternative,
                             mu = mu,
                             method_pval = method_pval,
                             ...)

    ## analysis configuration
    rst_conf <- list(Method       = dta_psrst$Method,
                     Outcome_type = dta_psrst$Outcome_type,
                     Study_type   = ifelse(is_rct, "RCT", "single-arm"))
    rst_conf$CI <- dta_psrst$CI[c("Method_ci",
                                  "Conf_type",
                                  "Conf_int",
                                  "Conf_stderr")]
    rst_conf$INFER <- dta_psrst$INFER[c("Method_infer",
                                        "Alternative",
                                        "Mu",
                                        "Method_pval")]
    ## summary observed
    rst_obs <- dta_psrst$Observed
    rst_obs$Group <- factor(rst_obs$Group,
                            levels = c(0, 1),
                            labels = c("RWD", "Cur"))
    if (is_rct) {
        rst_obs$Arm <- factor(rst_obs$Arm,
                              levels = c(0, 1),
                              labels = c("ctl", "trt"))
    } else {
        rst_obs$Arm <- NULL
    }

    ## summary estimation
    if (dta_psrst$Method %in% c("ps_km", "ps_lrk", "ps_rmst")) {
        ## i.e., ps_km, ps_lrk, and ps_rmst
        id_s <- factor(c(dta_psrst[[type]]$Stratum_Estimate$Stratum,
                         rep(0, nrow(dta_psrst[[type]]$Overall_Estimate))),
                       levels = c(1:nrow(dta_psrst$Borrow), 0),
                       labels = c(dta_psrst$Borrow$Stratum, "Overall"))
    } else {
        ## i.e., ps_pp and ps_cl
        id_s <- factor(c(1:nrow(dta_psrst$Borrow), 0),
                       levels = c(1:nrow(dta_psrst$Borrow), 0),
                       labels = c(dta_psrst$Borrow$Stratum, "Overall"))
    }
    dtype <- rbind(dta_psrst[[type]]$Stratum_Estimate[, col_est],
                   dta_psrst[[type]]$Overall_Estimate[, col_est])
    rst_est <- data.frame(Stratum = id_s)
    rst_est <- cbind(rst_est, dtype)
    rst_est_trt <- NULL
    rst_est_ctl <- NULL

    if (is_rct) {
        if (exists("Treatment", dta_psrst) && !is.null(dta_psrst$Treatment)) {
            dtype <- rbind(dta_psrst$Treatment$Stratum_Estimate[, col_est],
                           dta_psrst$Treatment$Overall_Estimate[, col_est])
            rst_est_trt <- data.frame(Stratum = id_s)
            rst_est_trt <- cbind(rst_est_trt, dtype)
        }

        if (exists("Control", dta_psrst) && !is.null(dta_psrst$Control)) {
            dtype <- rbind(dta_psrst$Control$Stratum_Estimate[, col_est],
                           dta_psrst$Control$Overall_Estimate[, col_est])
            rst_est_ctl <- data.frame(Stratum = id_s)
            rst_est_ctl <- cbind(rst_est_ctl, dtype)
        }
    }

    ## summary CI results
    dtype <- rbind(dta_psrst$CI[[type]]$Stratum_Estimate,
                   dta_psrst$CI[[type]]$Overall_Estimate)
    rst_est <- cbind(rst_est, dtype)

    if (is_rct) {
        if (exists("Treatment", dta_psrst) && !is.null(dta_psrst$Treatment)) {
            dtype <- rbind(dta_psrst$CI$Treatment$Stratum_Estimate,
                           dta_psrst$CI$Treatment$Overall_Estimate)
            rst_est_trt <- cbind(rst_est_trt, dtype)
        }

        if (exists("Control", dta_psrst) && !is.null(dta_psrst$Control)) {
            dtype <- rbind(dta_psrst$CI$Control$Stratum_Estimate,
                           dta_psrst$CI$Control$Overall_Estimate)
            rst_est_ctl <- cbind(rst_est_ctl, dtype)
        }
    }

    ## summary INFR results
    dtype <- rbind(dta_psrst$INFER[[type]]$Stratum_InferProb,
                   dta_psrst$INFER[[type]]$Overall_InferProb)
    if (dta_psrst$Method == "ps_pp") {
        colnames(dtype) <- "PostPr"
    } else {
        colnames(dtype) <- "p.value"
    }
    rst_est <- cbind(rst_est, dtype)

    ## return
    rst <- dta_psrst
    rst$OUTANA <- list(Analysis_Setup   = rst_conf,
                       Observed_Summary = rst_obs,
                       Analysis_Summary = rst_est)

    if (is_rct) {
        rst$OUTANA$RCT_Summary <- list(Treatment = rst_est_trt,
                                       Control   = rst_est_ctl)
    }

    class(rst) <- get_rwe_class("OUTANA")
    rst
}


#' @title Print outcome analysis results
#'
#' @description Print detailed information of outcome analysis results
#'
#' @param x A list of class \code{PSRWE_RST_OUTANA} that is generated using the
#'     \code{\link{psrwe_outana}} function.
#' @param show_details Print out more observed summary
#' @param show_rct Print out more analysis summary for RCT arms
#' @param show_pred_tps Specified time points to be shown
#' @param ... Additional parameters
#'
#'
#' @method print PSRWE_RST_OUTANA
#'
#' @return None (invisible \code{NULL})
#'
#' @export
#'
print.PSRWE_RST_OUTANA <- function(x,
                                   show_details = FALSE,
                                   show_rct = FALSE,
                                   show_pred_tps = NULL,
                                   ...) {

    ## check
    is_km <- x$Method %in% get_rwe_class("ANAMETHOD_KM")
    if (is_km) {
        if (is.null(show_pred_tps)) {
            pred_tps <- x$pred_tp
        } else {
            pred_tps <- show_pred_tps
        }
        x <- get_psrst_km_subset(x, pred_tps)
    }

    ## Detatch outana
    x_outana <- x$OUTANA

    ## Print study design
    cat(paste("- Method: ", x_outana$Analysis_Setup$Method,
              ", Outcome Type: ", x_outana$Analysis_Setup$Outcome_type,
              ", Study Type: ", x_outana$Analysis_Setup$Study_type,
              sep = ""))

    if (exists("pred_tp", x)) {
        cat(paste("\n- Predict Time Point: ",
                  paste(x$pred_tp, collapse = " "),
                  sep = ""))
    }

    if (exists("stderr_method", x)) {
        cat(paste("\n- StdErr Method: ",
                  x$stderr_method,
                  sep = ""))
    }

    cat("\n")

    if (exists("CI", x_outana$Analysis_Setup)) {
        ci <- x_outana$Analysis_Setup$CI
        cat(paste("- Interval Method: ", ci$Method_ci,
		  ", Level: ", ci$Conf_int,
                  sep = ""))

        if (!is.na(ci$Conf_type)) {
            cat(paste(", Type: ", ci$Conf_type,
                      sep = ""))
        }

        if (!is.na(ci$Conf_stderr)) {
            cat(paste(", Stderr: ", ci$Conf_stderr,
                      sep = ""))
        }


        cat("\n")
    }

    if (exists("INFER", x_outana$Analysis_Setup)) {
        infer <- x_outana$Analysis_Setup$INFER
        cat(paste("- Test Method: ", infer$Method_infer, sep = ""))

        if (infer$Method_infer == "p_value") {
            cat(paste(", Method pval: ", infer$Method_pval,
                      sep = ""))
        }

        cat("\n")

        if (x_outana$Analysis_Setup$Study_type == "RCT") {
            if (x_outana$Analysis_Setup$Method %in%
                c("ps_pp", "ps_cl", "ps_km")) {
                poi <- "theta_trt-theta_ctl"
            } else if(x_outana$Analysis_Setup$Method == "ps_lrk") {
                poi <- "sum[d_trt-E(d_trt)]"
            } else if(x_outana$Analysis_Setup$Method == "ps_rmst") {
                poi <- "auc(S_trt)-auc(S_ctl)"
            } else {
              stop("The method is not implemented.")
            }
        } else {
            poi <- "theta"
        }

        if (infer$Alternative == "less") {
            sign_null <- ">="
            sign_alter <- "<"
        } else if (infer$Alternative == "greater") {
            sign_null <- "<="
            sign_alter <- ">"
        } else {
            sign_null <- "=="
            sign_alter <- "!="
        }

        cat(paste("  H0: ", poi, " ", sign_null,
                  " ", sprintf("%5.3f", infer$Mu),
		  " vs. ", 
                  "Ha: ", poi, " ", sign_alter,
                  " ", sprintf("%5.3f", infer$Mu),
		  "\n", sep = ""))
    }

    ## Print summary statistics
    # cat("- Observed Data Summary:\n")
    # if (show_details) {
    #     print(x_outana$Observed_Summary, row.names = FALSE)
    # } else {
    #     print(x_outana$Observed_Summary[x_outana$Observed_Summary$Stratum ==
    #                                     "Overall",],
    #           row.names = FALSE)
    # }

    ## Print outcome analyses
    cat("- Analysis Results:\n")
    if (show_details) {
        print(x_outana$Analysis_Summary, row.names = FALSE)
    } else {
        print(x_outana$Analysis_Summary[x_outana$Analysis_Summary$Stratum ==
                                        "Overall",],
              row.names = FALSE)
    }

    if (exists("RCT_Summary", x_outana) && show_rct) {
        cat("- RCT Treatment Arm:\n")
        if (!is.null(x_outana$RCT_Summary$Treatment)) {
            if (show_details) {
                print(x_outana$RCT_Summary$Treatment, row.names = FALSE)
            } else {
                print(x_outana$RCT_Summary$Treatment[x_outana$RCT_Summary$Treatment$Stratum ==
                                                     "Overall",],
                      row.names = FALSE)
            }
        } else {
            print(NA)
        }

        cat("- RCT Control Arm:\n")
        if (!is.null(x_outana$RCT_Summary$Control)) {
            if (show_details) {
                print(x_outana$RCT_Summary$Control, row.names = FALSE)
            } else {
                print(x_outana$RCT_Summary$Control[x_outana$RCT_Summary$Control$Stratum ==
                                                   "Overall",],
                      row.names = FALSE)
            }
        } else {
            print(NA)
        }
    }

    return(invisible())
}


#' @title Summary outcome analysis results
#'
#' @description Summary information of outcome analysis results
#'
#' @param object A list of class \code{PSRWE_RST_OUTANA} that is generated using the
#'     \code{\link{psrwe_outana}} function.
#' @param pred_tps Specified time points
#' @param ... Additional parameters
#'
#'
#' @method summary PSRWE_RST_OUTANA
#'
#' @return A list of class \code{PSRWE_RST_OUTANA} with additiona information
#'
#' @export
#'
summary.PSRWE_RST_OUTANA <- function(object,
                                     pred_tps = NULL,
                                     ...) {
    ## check
    is_km <- object$Method %in% get_rwe_class("ANAMETHOD_KM")
    if (is_km) {
        if (is.null(pred_tps)) {
            pred_tps <- object$pred_tp
        }
        object <- get_psrst_km_subset(object, pred_tps)
    }

    return(object)
}


#' @title Obtain subset of psrst for given time point (pred_tp)
#'
#' @noRd
get_psrst_km_subset <- function(dta_psrst, pred_tps = NULL) {
    ## check components
    is_rct <- dta_psrst$is_rct
    if (is_rct) {
        types_est <- c("Control", "Treatment", "Effect")
    } else {
        types_est <- c("Control")
    }

    ## find closest time point (infimum)
    org_pred_tps <- dta_psrst$pred_tp
    if (!is.null(pred_tps)) {
        unique_tps <- sort(unique(dta_psrst$Observed$T))
        unique_pred_tps <- sort(unique(pred_tps))

        stopifnot(min(unique_pred_tps) >= min(unique_tps))

        org_pred_tps <- NULL
        for (i_tp in unique_pred_tps) {
            id <- which(unique_tps <= i_tp)
            org_pred_tps <- c(org_pred_tps, unique_tps[id[length(id)]])
        }

        time_table <- data.frame(new_T = unique_pred_tps, org_T = org_pred_tps)

        dta_psrst$pred_tp <- unique_pred_tps
    }

    ## index Observed and Control
    id_Observed_T <- dta_psrst$Observed$T %in% org_pred_tps
    id_Stratum_T <- dta_psrst$Control$Overall_Estimate$T %in% org_pred_tps
    id_Overall_T <- dta_psrst$Control$Overall_Estimate$T %in% org_pred_tps

    ## Subset df by id_T and replace T with new time points in t_tbl
    subset_replace <- function(df, id_T, t_tbl) {
        df <- df[id_T,]
        if ("T" %in% names(df)) {
            tmp_T <- df$T
            for (i_t in 1:nrow(t_tbl)) {
              tmp_T[df$T == t_tbl$org_T[i_t]] <- t_tbl$new_T[i_t]
            }
            df$T <- tmp_T
        }
        return(df)
    }

    ## subset Observed
    dta_psrst$Observed <- subset_replace(dta_psrst$Observed,
                                         id_Observed_T,
                                         time_table)

    ## subset estimates
    for (i_type in types_est) {
        if (!is.null(dta_psrst[[i_type]])) {
            dta_psrst[[i_type]]$Stratum_Estimate <-
                subset_replace(dta_psrst[[i_type]]$Stratum_Estimate,
                               id_Stratum_T,
                               time_table)

            dta_psrst[[i_type]]$Overall_Estimate <-
                subset_replace(dta_psrst[[i_type]]$Overall_Estimate,
                               id_Overall_T,
                               time_table)
        }
    }

    ## subset CI
    for (i_type in types_est) {
        if (!is.null(dta_psrst$CI[[i_type]])) {
            dta_psrst$CI[[i_type]]$Stratum_Estimate <-
                subset_replace(dta_psrst$CI[[i_type]]$Stratum_Estimate,
                               id_Stratum_T,
                               time_table)

            dta_psrst$CI[[i_type]]$Overall_Estimate <-
                subset_replace(dta_psrst$CI[[i_type]]$Overall_Estimate,
                               id_Overall_T,
                               time_table)
        }
    }

    ## subset INFR
    for (i_type in types_est) {
        if (!is.null(dta_psrst$INFR[[i_type]])) {
            dta_psrst$INFR[[i_type]]$Stratum_Estimate <-
                subset_replace(dta_psrst$INFR[[i_type]]$Stratum_Estimate,
                               id_Stratum_T,
                               time_table)

            dta_psrst$INFR[[i_type]]$Overall_Estimate <-
                subset_replace(dta_psrst$INFR[[i_type]]$Overall_Estimate,
                               id_Overall_T,
                               time_table)
        }
    }

    ## subset OUTANA
    if (exists("OUTANA", dta_psrst)) {
        id_T <- dta_psrst$OUTANA$Observed_Summary$T %in% org_pred_tps
        dta_psrst$OUTANA$Observed_Summary <-
            subset_replace(dta_psrst$OUTANA$Observed_Summary,
                           id_T,
                           time_table)

        id_T <- dta_psrst$OUTANA$Analysis_Summary$T %in% org_pred_tps
        dta_psrst$OUTANA$Analysis_Summary <-
            subset_replace(dta_psrst$OUTANA$Analysis_Summary,
                           id_T,
                           time_table)

        if (is_rct) {
            id_T <- dta_psrst$OUTANA$RCT_Summary$Treatment$T %in% org_pred_tps
            dta_psrst$OUTANA$RCT_Summary$Treatment <-
                subset_replace(dta_psrst$OUTANA$RCT_Summary$Treatment,
                               id_T,
                               time_table)

            id_T <- dta_psrst$OUTANA$RCT_Summary$Control$T %in% org_pred_tps
            dta_psrst$OUTANA$RCT_Summary$Control <-
                subset_replace(dta_psrst$OUTANA$RCT_Summary$Control,
                               id_T,
                               time_table)
        }
    }

    return(dta_psrst)
}

