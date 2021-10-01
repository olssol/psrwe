#' @title Outcome Analysis for PS-Integrated Estimation
#'
#' Report outcome analysis for the PS-integrated approach.
#'
#' @param dta_psrst a returned object with class \code{PSRWE_EST}
#' @param method_ci a method name for confidence interval (default Wald)
#' @param conf_type a type name of transformation for the confidence interal
#'        of PSKM approach (default log_log)
#' @param conf_int a two-sided level of confidence/credible limits
#'        (default 0.95)
#' @param alternative a character string for the alternative hypothesis that
#'        must be one of \code{"less"} (default) or \code{"greater"}
#' @param mu a number indicating the true value of the parameter of interest
#'        (or the difference in means for two arms)
#' @param ... other options
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
                         alternative = c("less", "greater"),
                         mu = 0,
                         ...) {
    ## check
    stopifnot(inherits(dta_psrst,
                       what = get_rwe_class("ANARST")))

    stopifnot(dta_psrst$Method %in% c("ps_pp", "ps_cl", "ps_km"))

    ## check components
    outcome_type <- dta_psrst$Outcome_type
    is_rct <- dta_psrst$is_rct
    if (is_rct) {
        type <- "Effect"
    } else {
        type <- "Control"
    }
    is_km <- dta_psrst$Method == "ps_km"

    ## add ci and infer
    dta_psrst <- psrwe_ci(dta_psrst,
                          method_ci = method_ci[1],
                          conf_type = conf_type[1],
                          conf_int = conf_int[1],
                          ...)
    dta_psrst <- psrwe_infer(dta_psrst,
                             alternative = alternative[1],
                             mu = mu[1])

    ## analysis configuration
    rst_conf <- list(Method       = dta_psrst$Method,
                     Outcome_type = dta_psrst$Outcome_type,
                     Study_type   = ifelse(is_rct, "RCT", "single-arm"))
    if (is_km) {
        rst_conf$pred_tp <- dta_psrst$pred_tp
    }
    rst_conf$CI <- dta_psrst$CI[c("Method_ci",
                                  "Conf_type",
                                  "Conf_int")]
    rst_conf$INFER <- dta_psrst$INFER[c("Method_infer",
                                        "Alternative",
                                        "Mu")]

    ## summary observed
    dtype <- dta_psrst$Observed
    if (is_km) {
        dtype <- data.frame(dtype) %>%
            filter(dta_psrst$pred_tp == T)
        dtype <- data.frame(dtype[, colnames(dtype) != "T"])
    }
    rst_obs <- dtype
    if (!is_km) {
        colnames(rst_obs)[colnames(rst_obs) == "StdErr"] <- "SD"
    }

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
    col_est <- c("Mean", "StdErr")
    if (is_km) {
        col_est <- c(col_est, "T")
    }

    dtype <- rbind(dta_psrst[[type]]$Stratum_Estimate[, col_est],
                   dta_psrst[[type]]$Overall_Estimate[, col_est])
    if (is_km) {
        dtype <- data.frame(dtype) %>% filter(dta_psrst$pred_tp == T)
        dtype <- data.frame(dtype[, colnames(dtype) != "T"])
    }
    rst_est <- data.frame(Stratum = c(dta_psrst$Borrow$Stratum,
                                      "Overall"))
    rst_est <- cbind(rst_est, dtype)

    if (is_rct) {
        dtype <- rbind(dta_psrst$Treatment$Stratum_Estimate[, col_est],
                       dta_psrst$Treatment$Overall_Estimate[, col_est])
        if (is_km) {
            dtype <- data.frame(dtype) %>%
                filter(dta_psrst$pred_tp == T)
            dtype <- data.frame(dtype[, colnames(dtype) != "T"])
        }
        rst_est_trt <- data.frame(Stratum = c(dta_psrst$Borrow$Stratum,
                                          "Overall"))
        rst_est_trt <- cbind(rst_est_trt, dtype)

        dtype <- rbind(dta_psrst$Control$Stratum_Estimate[, col_est],
                       dta_psrst$Control$Overall_Estimate[, col_est])
        if (is_km) {
            dtype <- data.frame(dtype) %>%
                filter(dta_psrst$pred_tp == T)
            dtype <- data.frame(dtype[, colnames(dtype) != "T"])
        }
        rst_est_ctl <- data.frame(Stratum = c(dta_psrst$Borrow$Stratum,
                                          "Overall"))
        rst_est_ctl <- cbind(rst_est_ctl, dtype)
    }

    ## summary CI results
    dtype <- rbind(dta_psrst$CI[[type]]$Stratum_Estimate,
                   dta_psrst$CI[[type]]$Overall_Estimate)
    if (is_km) {
        dtype <- cbind(dtype,
                       T = c(dta_psrst[[type]]$Stratum_Estimate$T,
                             dta_psrst[[type]]$Overall_Estimate$T))
        dtype <- data.frame(dtype) %>%
            filter(dta_psrst$pred_tp == T)
        dtype <- data.frame(dtype[, colnames(dtype) != "T"])
    }
    rst_est <- cbind(rst_est, dtype)

    if (is_rct) {
        dtype <- rbind(dta_psrst$CI$Treatment$Stratum_Estimate,
                       dta_psrst$CI$Treatment$Overall_Estimate)
        if (is_km) {
            dtype <- cbind(dtype,
                           T = c(dta_psrst$Treatment$Stratum_Estimate$T,
                                 dta_psrst$Treatment$Overall_Estimate$T))
            dtype <- data.frame(dtype) %>%
                filter(dta_psrst$pred_tp == T)
            dtype <- data.frame(dtype[, colnames(dtype) != "T"])
        }
        rst_est_trt <- cbind(rst_est_trt, dtype)

        dtype <- rbind(dta_psrst$CI$Control$Stratum_Estimate,
                       dta_psrst$CI$Control$Overall_Estimate)
        if (is_km) {
            dtype <- cbind(dtype,
                           T = c(dta_psrst$Control$Stratum_Estimate$T,
                                 dta_psrst$Control$Overall_Estimate$T))
            dtype <- data.frame(dtype) %>%
                filter(dta_psrst$pred_tp == T)
            dtype <- data.frame(dtype[, colnames(dtype) != "T"])
        }
        rst_est_ctl <- cbind(rst_est_ctl, dtype)
    }

    ## summary INFR results
    dtype <- rbind(dta_psrst$INFER[[type]]$Stratum_InferProb,
                   dta_psrst$INFER[[type]]$Overall_InferProb)
    if (is_km) {
        dtype <- cbind(dtype,
                       T = c(dta_psrst[[type]]$Stratum_Estimate$T,
                             dta_psrst[[type]]$Overall_Estimate$T))
        dtype <- data.frame(dtype) %>%
            filter(dta_psrst$pred_tp == T)
        dtype <- data.frame(dtype[, colnames(dtype) != "T"])
    }

    if (dta_psrst$Method == "ps_pp") {
        colnames(dtype) <- "PostPr"
    } else {
        colnames(dtype) <- "p-value"
    }
    rst_est <- cbind(rst_est, dtype)

    ## return
    rst <- dta_psrst
    rst$OUTANA <- list(Analysis_Setup  = rst_conf,
                       Observed_Summary  = rst_obs,
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
#' @description Print summary information of outcome analysis results
#'
#' @param x A list of class \code{PSRWE_RST_OUTANA} that is generated using the
#'     \code{\link{psrwe_outana}} function.
#' @param show_details print out more observed summary
#' @param show_rct print out more analysis summary for RCT arms
#' @param ... Additional parameters
#'
#'
#' @method print PSRWE_RST_OUTANA
#'
#'
#' @export
#'
print.PSRWE_RST_OUTANA <- function(x,
                                   show_details = FALSE,
                                   show_rct = FALSE,
                                   ...) {
    ## Detatch outana
    x_outana <- x$OUTANA

    ## Print study design
    cat(paste("- Method: ", x_outana$Analysis_Setup$Method,
              ", Outcome Type: ", x_outana$Analysis_Setup$Outcome_type,
              ", Study Type: ", x_outana$Analysis_Setup$Study_type,
              sep = ""))

    if (exists("pred_tp", x_outana$Analysis_Setup)) {
        cat(paste("\n- Predict Time Point: ",
                  x_outana$Analysis_Setup$pred_tp,
                  sep = ""))
    }

    cat("\n")

    if (exists("CI", x_outana$Analysis_Setup)) {
        ci <- x_outana$Analysis_Setup$CI
        cat(paste("- Confidence Interval: ", ci$Method_ci,
		  ", Level: ", ci$Conf_int,
                  sep = ""))

        if (!is.na(ci$Conf_type)) {
            cat(paste(", Type: ", ci$Conf_type,
                      sep = ""))
        }

        cat("\n")
    }

    if (exists("INFER", x_outana$Analysis_Setup)) {
        infer <- x_outana$Analysis_Setup$INFER
        cat(paste("- Test Method: ", infer$Method_infer,
		  "\n", sep = ""))

        if (x_outana$Analysis_Setup$Study_type == "RCT") {
            poi <- "theta_trt-theta_ctl"
        } else {
            poi <- "theta"
        }

        cat(paste("  H0: ", poi, " ",
		  ifelse(infer$Alternative == "less", ">=", "<="),
                  " ", sprintf("%5.3f", infer$Mu),
		  " vs. ", 
                  "Ha: ", poi, " ",
		  ifelse(infer$Alternative == "less", "<", ">"),
                  " ", sprintf("%5.3f", infer$Mu),
		  "\n", sep = ""))
    }

    ## Print summary statistics
    cat("- Observed Data Summary:\n")
    if (show_details) {
        print(x_outana$Observed_Summary)
    } else {
        print(x_outana$Observed_Summary[x_outana$Observed_Summary$Stratum ==
                                        "Overall",])
    }

    ## Print outcome analyses
    cat("- Analysis Results:\n")
    print(x_outana$Analysis_Summary)

    if (exists("RCT_Summary", x_outana) && show_rct) {
        cat("- RCT Treatment Arm:\n")
        print(x_outana$RCT_Summary$Treatment)

        cat("- RCT Control Arm:\n")
        print(x_outana$RCT_Summary$Control)
    }
    invisible()
}

