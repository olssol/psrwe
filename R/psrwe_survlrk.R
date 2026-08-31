#' PS-Integrated Log-Rank Test For Comparing Time-to-event Outcomes
#'
#' Log-rank test evaluates a two-arm RCT for up to a given time point.
#' Variance can be estimated by Jackknife methods.
#' Apply to the case when there is only one external data source and
#' two-arm RCT.
#'
#' @inheritParams psrwe_survkm
#'
#' @param pred_tp A numeric value corresponding to the time of interest
#'                (e.g., 365 days or 1 year)
#' @param v_time Column name corresponding to event time
#' @param v_event Column name corresponding to event status
#' @param stderr_method Method for computing StdErr (see Details)
#' @param n_bootstrap Number of bootstrap samples (for bootstrap stderr)
#' @param ... Additional Parameters
#'
#' @details \code{stderr_method} includes \code{naive} as default which
#'     mostly follows the Greenwood formula,
#'     \code{jk} using the Jackknife method within each stratum,
#'     \code{sjk} using simple Jackknife method for combined estimates
#'     such as point estimates in single-arm or treatment effects in RCT, or
#'     \code{cjk} for complex Jackknife method including refitting PS model,
#'     matching, trimming, calculating borrowing parameters, and
#'     combining overall estimates.
#'     Note that \code{sjk} may take a while longer to finish and
#'     \code{cjk} will take even longer to finish.
#'     The \code{sbs} and \code{cbs} is for simple and complex Bootstrap
#'     methods.
#'
#' @return A data frame with class name \code{PSRWE_RST_TESTANA}.
#'     It contains the test statistics of each stratum as well as the
#'     Jackknife estimation. The results can be further
#'     summarized by its S3 method \code{summary}.
#'     The results can also be analyzed by \code{psrwe_outana} for outcome
#'     analysis and inference.
#'
#'
#' @examples
#' data(ex_dta_rct)
#' dta_ps_rct <- psrwe_est(ex_dta_rct,
#'                         v_covs = paste("V", 1:7, sep = ""),
#'                         v_grp = "Group", cur_grp_level = "current",
#'                         v_arm = "Arm", ctl_arm_level = "control",
#'                         ps_method = "logistic", nstrata = 5,
#'                         stra_ctl_only = FALSE)
#' ps_bor_rct <- psrwe_borrow(dta_ps_rct, total_borrow = 30)
#' rst_lrk <- psrwe_survlrk(ps_bor_rct,
#'                          pred_tp = 365, 
#'                          v_time = "Y_Surv",
#'                          v_event = "Status")
#' rst_lrk
#'
#' @export
#'
psrwe_survlrk <- function(dta_psbor, pred_tp,
                          v_time        = "time",
                          v_event       = "event",
                          stderr_method = c("naive", "jk", "sjk", "cjk",
                                            "sbs", "cbs", "none"), 
                          n_bootstrap = 200,
                          ...) {

    ## check
    stopifnot(dta_psbor$is_rct)

    stopifnot(inherits(dta_psbor,
                       what = get_rwe_class("PSDIST")))

    stopifnot(all(c(v_event, v_time) %in%
                  colnames(dta_psbor$data)))

    stopifnot(is.numeric(pred_tp) & 1 == length(pred_tp))

    stderr_method <- match.arg(stderr_method)

    ## all time points
    data    <- dta_psbor$data
    obs_tps <- data[which(1 == data[[v_event]]), v_time]
    all_tps <- sort(unique(c(pred_tp, obs_tps)))

    ## observed
    rst_obs <- get_km_observed(data, v_time, v_event, all_tps)

    ## call estimation
    if (stderr_method %in% c("naive", "jk", "none")) {
        rst <- get_ps_lrk_rmst(dta_psbor,
                               v_event = v_event, v_time = v_time,
                               f_stratum = get_surv_stratum_lrk,
                               pred_tps = all_tps,
                               stderr_method = stderr_method,
                               ...)
    } else if (stderr_method %in% c("sjk")) {
        rst <- get_ps_lrk_rmst_sjk(dta_psbor,
                                   v_event = v_event, v_time = v_time,
                                   f_stratum = get_surv_stratum_lrk,
                                   pred_tps = all_tps,
                                   stderr_method = "none",
                                   ...)
    } else if (stderr_method %in% c("cjk")) {
        rst <- get_ps_lrk_rmst_cjk(dta_psbor,
                                   v_event = v_event, v_time = v_time,
                                   f_stratum = get_surv_stratum_lrk,
                                   pred_tp = all_tps,
                                   stderr_method = "none",
                                   ...)
    } else if (stderr_method %in% c("sbs")) {
        rst <- get_ps_lrk_rmst_sbs(dta_psbor,
                                   v_event = v_event, v_time = v_time,
                                   f_stratum = get_surv_stratum_lrk,
                                   pred_tps = all_tps,
                                   stderr_method = "none",
                                   n_bootstrap = n_bootstrap,
                                   ...)
    } else if (stderr_method %in% c("cbs")) {
        rst <- get_ps_lrk_rmst_cbs(dta_psbor,
                                   v_event = v_event, v_time = v_time,
                                   f_stratum = get_surv_stratum_lrk,
                                   pred_tp = all_tps,
                                   stderr_method = "none",
                                   n_bootstrap = n_bootstrap,
                                   ...)
    } else {
        stop("stderr_errmethod is not implemented.")
    }

    ## return
    rst$Observed <- rst_obs
    rst$pred_tp  <- pred_tp
    rst$stderr_method <- stderr_method
    rst$Method   <- "ps_lrk"
    rst$Outcome_type <- "tte"
    class(rst)   <- get_rwe_class("ANARST")
    return(rst)
}

#' Get log-rank estimation for each stratum
#'
#'
#' @noRd
#'
get_surv_stratum_lrk <- function(d1, d0 = NULL, d1t, n_borrow = 0, pred_tps,
                                 stderr_method = "jk", ...) {

    ## treatment or control only
    dta_cur <- d1
    dta_ext <- d0
    dta_cur_trt <- d1t
    ns1     <- nrow(dta_cur)
    ns1_trt <- nrow(dta_cur_trt)

    if (is.null(d0)) {
        ns0 <- 0
    } else {
        ns0 <- nrow(dta_ext)
    }

    ## overall estimate
    overall  <- rwe_lrk(dta_cur, dta_ext, dta_cur_trt, n_borrow, pred_tps,
                        stderr_method)

    ## jackknife stderr
    if (stderr_method == "jk") {
        overall_theta <- overall[, 1, drop = TRUE]

        jk_theta      <- rep(0, length(overall_theta))
        for (j in seq_len(ns1)) {
            cur_jk   <- rwe_lrk(dta_cur[-j, ], dta_ext, dta_cur_trt, n_borrow,
                                pred_tps, stderr_method)
            jk_theta <- jk_theta + (cur_jk[, 1] - overall_theta)^2
        }

        for (j in seq_len(ns1_trt)) {
            cur_jk   <- rwe_lrk(dta_cur, dta_ext, dta_cur_trt[-j, ], n_borrow,
                                pred_tps, stderr_method)
            jk_theta <- jk_theta + (cur_jk[, 1] - overall_theta)^2
        }

        if (ns0 > 0) {
            for (j in seq_len(ns0)) {
                ext_jk   <- rwe_lrk(dta_cur, dta_ext[-j, ], dta_cur_trt,
                                    n_borrow, pred_tps, stderr_method)
                jk_theta <- jk_theta + (ext_jk[, 1] - overall_theta)^2
            }
        }

        ## summary
        nc_jk <- (ns1 + ns0 + ns1_trt - 1) / (ns1 + ns0 + ns1_trt)
        stderr_theta <- sqrt(jk_theta * nc_jk)

        overall[, 2] <- stderr_theta
    }

    return(overall)
}

#' Log-rank Estimation
#'
#' Estimate log-rank estimates for a single PS
#' stratum
#'
#'
#' @param dta_cur Matrix of time and event from a PS stratum in the current
#'                study (control arm only)
#' @param dta_ext Matrix of time and event from a PS stratum in the external
#'                data source (control arm only)
#' @param dta_cur_trt Matrix of time and event from a PS stratum in current
#'                    study (treatment arm only)
#' @param n_borrow Number of subjects to be borrowed
#' @param pred_tps All time points of events (unique and sorted)
#' @param stderr_method Method for computing StdErr (available for naive only)
#'
#' @return Estimation of log-rank estimates at time \code{pred_tps}
#'
#'
#' @export
#'
rwe_lrk <- function(dta_cur, dta_ext, dta_cur_trt, n_borrow = 0,
                    pred_tps = NULL, stderr_method = "naive") {

    ## current control and external control if available
    cur_data    <- dta_cur
    ns1         <- nrow(dta_cur)
    cur_weights <- rep(1, ns1)

    if (n_borrow > 0) {
        ns0         <- nrow(dta_ext)
        cur_data    <- rbind(cur_data, dta_ext)
        cur_weights <- c(cur_weights,
                         rep(n_borrow / ns0, ns0))
    }

    ## KM with stratum weights
    colnames(cur_data) <- c("time", "event")
    cur_data <- data.frame(cur_data)
    cur_surv <- survfit(Surv(time, event) ~ 1,
                        data      = cur_data,
                        weights   = cur_weights,
                        conf.type = "none")

    ## trt arm
    cur_data_trt    <- dta_cur_trt
    ns1_trt         <- nrow(dta_cur_trt)
    cur_weights_trt <- rep(1, ns1_trt)

    ## KM with stratum weights for trt arm
    colnames(cur_data_trt) <- c("time", "event")
    cur_data_trt <- data.frame(cur_data_trt)
    cur_surv_trt <- survfit(Surv(time, event) ~ 1,
                            data      = cur_data_trt,
                            weights   = cur_weights_trt,
                            conf.type = "none")

    ## summary.survfit() need to be extend to longer time points
    ## Last values will be carried over for predictions
    if (is.null(pred_tps)) {
        pred_tps <- sort(unique(c(cur_surv$time[cur_surv$n.event > 0],
                                  cur_surv_trt$time[cur_surv_trt$n.event > 0])))
    }

    ## summary.survfit() only reports information for specified time points
    ## so, pred_tps should be inputted will all time points before tau.
    rst <- summary(cur_surv, time = pred_tps, extend = TRUE)
    rst_trt <- summary(cur_surv_trt, time = pred_tps, extend = TRUE)

    ## info needed for LRK
    n_risk_trt <- rst_trt$n.risk
    n_risk_ctl <- rst$n.risk
    n_event_trt <- rst_trt$n.event
    n_event_ctl <- rst$n.event

    n_risk <- n_risk_trt + n_risk_ctl
    n_event <- n_event_trt + n_event_ctl
    p_event <- ifelse(n_risk == 0, 0, n_event / n_risk)
    E_1_j <- n_risk_trt * p_event

    ## LRK main statistic
    mean_d <- n_event_trt - E_1_j
    mean_d <- cumsum(mean_d)

    ## for LRK naive stderr
    if (stderr_method == "naive") {
        stderr_d <- ifelse(n_risk <= 1, 0,
                           E_1_j * (1 - p_event) *
                           n_risk_ctl / (n_risk - 1))
        stderr_d <- sqrt(cumsum(stderr_d))
    } else {
        ## for none, jk, sjk, cjk, sbs, or cbs
        stderr_d <- rep(NA, length(mean_d))
    }

    ## combine log-rank estimates
    rst_lrk <- cbind(mean_d, stderr_d, pred_tps)

    colnames(rst_lrk) <- c("Mean", "StdErr", "T")
    return(rst_lrk)
}

