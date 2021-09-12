#' @title Outcome Analysis for PS-Integrated Estimation
#'
#' Report outcome analysis for the PS-integrated approach.
#'
#' @param dta_psrst a returned object with class \code{PS_RWE_EST}
#' @param ... other options
#'
#' @return A list with class name \code{PS_RWE_EST_OUTANA}.
#'
#' @details This function is mainly for summarizing and reporting the
#'     outcome analysis for the PS-integrated estimation.
#'     The input \code{dta_psrst} can be generated from the functions
#'     \code{\link{ps_rwe_powerp}}, \code{\link{ps_rwe_compl}},
#'     \code{\link{ps_rwe_survkm}}, \code{\link{ps_rwe_ci}}, and.
#'     \code{\link{ps_rwe_infer}}.
#'
#' @examples
#' data(ex_dta)
#' dta_ps <- ps_rwe_est(ex_dta,
#'        v_covs = paste("V", 1:7, sep = ""),
#'        v_grp = "Group",
#'        cur_grp_level = "current")
#' ps_borrow <- ps_rwe_borrow(total_borrow = 30, dta_ps)
#' ps_rst <- ps_rwe_compl(ps_borrow, v_outcome = "Y_Con")
#' rst <- ps_rwe_outana(ps_rst)
#' rst
#'
#' @export
#'
ps_rwe_outana <- function(dta_psrst) {
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

    ## unique names
    u_group <- sort(unique(dta_psrst$Observed$Group))
    u_arm <- sort(unique(dta_psrst$Observed$Arm))
    u_stratum <- sort(dta_psrst$Borrow$Stratum)

    ## study configuration
    rst_conf <- list(Method       = dta_psrst$Method,
                     Outcome_type = dta_psrst$Outcome_type)
    if (is_km) {
        rst_conf$pred_tp <- dta_psrst$pred_tp
    }
    if (is_ci) {
       rst_conf <- c(rst_conf,
                     dta_psrst$CI[c("Method_ci",
                                    "Conf_type",
                                    "Conf_int")])
    }
    if (is_infer) {
       rst_conf <- c(rst_conf,
                     dta_psrst$INFER[c("Method_infer",
                                       "Alternative",
                                       "Mu")])
    }

    ## summary observed
    dtype <- dta_psrst$Observed
    if (is_km) {
        dtype <- data.frame(dtype) %>%
            filter(dta_psrst$pred_tp == T)
    }

    id <- NULL
    for (g in u_group) {
        for (a in u_arm) {
            for (s in u_stratum) {
                id_tmp <- which(dtype$Group   == g &
                                dtype$Arm     == a &
                                dtype$Stratum == s)
                id <- c(id, id_tmp)
            }
        }
    }
    for (g in u_group) {
        for (a in u_arm) {
            id_tmp <- which(dtype$Group   == g &
                            dtype$Arm     == a &
                            dtype$Stratum == "Overall")
            id <- c(id, id_tmp)
        }
    }

    rst_obs <- dtype[id, c("Group", "Arm", "Stratum", "N", "Mean", "StdErr")]

    ## summary estimation
    dtype <- rbind(dta_psrst[[type]]$Stratum_Estimate,
                   dta_psrst[[type]]$Overall_Estimate)
    if (is_km) {
        dtype <- data.frame(dtype) %>%
            filter(dta_psrst$pred_tp == T)
    }

    rst_est <- data.frame(Stratum = c(dta_psrst$Borrow$Stratum,
                                      "Overall"),
                          Mean    = dtype[, "Mean"],
                          StdErr  = dtype[, "StdErr"])

    if (is_rct) {
        colnames(rst_est) <- c("Stratum", "Mean_Eff", "StdErr_Eff")

        dtype <- rbind(dta_psrst$Control$Stratum_Estimate,
                       dta_psrst$Control$Overall_Estimate)
        if (is_km) {
            dtype <- data.frame(dtype) %>%
                filter(dta_psrst$pred_tp == T)
        }
        rst_est_ctl <- data.frame(Mean_CTL   = dtype[, "Mean"],
                                  StdErr_CTL = dtype[, "StdErr"])

        dtype <- rbind(dta_psrst$Treatment$Stratum_Estimate,
                       dta_psrst$Treatment$Overall_Estimate)
        if (is_km) {
            dtype <- data.frame(dtype) %>%
                filter(dta_psrst$pred_tp == T)
        }
        rst_est_trt <- data.frame(Mean_TRT   = dtype[, "Mean"],
                                  StdErr_TRT = dtype[, "StdErr"])

        rst_est <- cbind(rst_est,
                         rst_est_trt,
			 rst_est_ctl)
    }

    ## summary CI results
    if (is_ci) {
    }

    ## summary INFR results
    if (is_infer) {
    }


    ## return
    rst <- c(rst,
             list(Observed_Summary       = rst_obs,
                  Estimates_Summary      = rst_est,
                  Confidence_Interval    = rst_ci,
                  Inferential_Statistics = rst_infer))
    class(rst) <- get_rwe_class("OUTANA")
    rst
}


#' @title Print outcome analysis results
#'
#' @description Print summary information of outcome analysis results
#'
#' @param x A list of class \code{PS_RWE_RST_OUTANA} that is generated using the
#'     \code{\link{ps_rwe_outana}} function.
#' @param ... Additional parameters
#'
#'
#' @method print PS_RWE_RST_OUTANA
#'
#'
#' @export
#'
print.PS_RWE_RST_OUTANA <- function(x, ...) {
    print(x)
    invisible()
}

