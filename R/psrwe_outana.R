#' @title Outcome Analysis for PS-Integrated Estimation
#'
#' Report outcome analysis for the PS-integrated approach.
#'
#' @param dta_psrst a returned object with class \code{PSRWE_EST}
#' @param ... other options
#'
#' @return A list with class name \code{PSRWE_EST_OUTANA}.
#'
#' @details This function is mainly for summarizing and reporting the
#'     outcome analysis for the PS-integrated estimation.
#'     The input \code{dta_psrst} can be generated from the functions
#'     \code{\link{psrwe_powerp}}, \code{\link{psrwe_compl}},
#'     \code{\link{psrwe_survkm}}, \code{\link{psrwe_ci}}, and.
#'     \code{\link{psrwe_infer}}.
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
psrwe_outana <- function(dta_psrst) {
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
    is_ci <- exists("CI", dta_psrst)
    is_infer <- exists("INFER", dta_psrst)

    ## analysis configuration
    rst_conf <- list(Method       = dta_psrst$Method,
                     Outcome_type = dta_psrst$Outcome_type)
    if (is_km) {
        rst_conf$pred_tp <- dta_psrst$pred_tp
    }
    if (is_ci) {
       rst_conf$CI <- dta_psrst$CI[c("Method_ci",
                                     "Conf_type",
                                     "Conf_int")]
    }
    if (is_infer) {
       rst_conf$INFER <- dta_psrst$INFER[c("Method_infer",
                                           "Alternative",
                                           "Mu")]
    }

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
    if (is_ci) {
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
    }

    ## summary INFR results
    if (is_infer) {
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
    }


    ## return
    rst <- list(Analysis_Setup  = rst_conf,
                Observed_Summary  = rst_obs,
                Analysis_Summary = rst_est)

    if (is_rct) {
        rst$RCT_Summary <- list(Treatment = rst_est_trt,
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
    cat(paste("- Method: ", x$Analysis_Setup$Method,
              ", Outcome Type: ", x$Analysis_Setup$Outcome_type,
              sep = ""))

    if (exists("pred_tp", x$Analysis_Setup)) {
        cat(paste(", Predict Time: ",
                  x$Analysis_Setup$pred_tp,
                  sep = ""))
    }

    cat("\n")

    if (exists("CI", x$Analysis_Setup)) {
        ci <- x$Analysis_Setup$CI
        cat(paste("- Confidence Interval: ", ci$Method_ci,
		  ", Level: ", ci$Conf_int,
                  sep = ""))

        if (!is.na(ci$Conf_type)) {
            cat(paste(", Type: ", ci$Conf_type,
                      sep = ""))
        }

        cat("\n")
    }

    if (exists("INFER", x$Analysis_Setup)) {
        infer <- x$Analysis_Setup$INFER
        cat(paste("- Test Method: ", infer$Method_infer,
		  "\n", sep = ""))

        if (is_rct) {
            poi <- "theta_trt-theta_ctl"
        } else {
            poi <- "theta"
        }

        cat(paste("- H0: ", poi, " ",
		  ifelse(infer$Alternative == "less", ">=", "<="),
                  " ", sprintf("%5.3f", infer$Mu),
		  " vs. ", 
                  "Ha: ", poi, " ",
		  ifelse(infer$Alternative == "less", "<", ">"),
                  " ", sprintf("%5.3f", infer$Mu),
		  "\n", sep = ""))
    }

    cat("- Observed Data Summary:\n")
    if (show_details) {
        print(x$Observed_Summary)
    } else {
        print(x$Observed_Summary[x$Observed_Summary$Stratum == "Overall",])
    }

    cat("- Analysis Results:\n")
    print(x$Analysis_Summary)

    if (exists("RCT_Summary", x) && show_rct) {
        cat("- RCT Treatment Arm:\n")
        print(x$RCT$Treatment)

        cat("- RCT Control Arm:\n")
        print(x$RCT$Control)
    }
    invisible()
}

