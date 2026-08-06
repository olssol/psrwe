#' PS-Integrated Restricted Mean Survival Time (RMST) Test For Comparing Time-to-event Outcomes
#'
#' RMST test evaluates a two-arm RCT for up to a given time point.
#' Variance can be estimated by the Jackknife method.
#' Apply to the case when there is only one external data source and
#' two-arm RCT.
#'
#' @inheritParams psrwe_survkm
#'
#' @param pred_tp A numeric value corresponding to the time of interest
#'                (e.g., 365 days or 1 year)
#' @param v_time Column name corresponding to event time
#' @param v_event Column name corresponding to event status
#' @param stderr_method Method for computing StdErr, see Details
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
#' rst_rmst <- psrwe_survrmst(ps_bor_rct,
#'                            pred_tp = 365,
#'                            v_time = "Y_Surv",
#'                            v_event = "Status")
#' rst_rmst
#'
#' @export
#'
psrwe_survrmst <- function(dta_psbor, pred_tp,
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
                               f_stratum = get_surv_stratum_rmst,
                               pred_tps = all_tps,
                               stderr_method = stderr_method,
                               ...)
    } else if (stderr_method %in% c("sjk")) {
        rst <- get_ps_lrk_rmst_sjk(dta_psbor,
                                   v_event = v_event, v_time = v_time,
                                   f_stratum = get_surv_stratum_rmst,
                                   pred_tps = all_tps,
                                   stderr_method = "none",
                                   ...)
    } else if (stderr_method %in% c("cjk")) {
        rst <- get_ps_lrk_rmst_cjk(dta_psbor,
                                   v_event = v_event, v_time = v_time,
                                   f_stratum = get_surv_stratum_rmst,
                                   pred_tp = all_tps,
                                   stderr_method = "none",
                                   ...)
    } else if (stderr_method %in% c("sbs")) {
        rst <- get_ps_lrk_rmst_sbs(dta_psbor,
                                   v_event = v_event, v_time = v_time,
                                   f_stratum = get_surv_stratum_rmst,
                                   pred_tps = all_tps,
                                   stderr_method = "none",
                                   n_bootstrap = n_bootstrap,
                                   ...)
    } else if (stderr_method %in% c("cbs")) {
        rst <- get_ps_lrk_rmst_cbs(dta_psbor,
                                   v_event = v_event, v_time = v_time,
                                   f_stratum = get_surv_stratum_rmst,
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
    rst$Method   <- "ps_rmst"
    rst$Outcome_type <- "tte"
    class(rst)   <- get_rwe_class("ANARST")
    return(rst)
}

#' Get RMST estimation for each stratum
#'
#'
#' @noRd
#'
get_surv_stratum_rmst <- function(d1, d0 = NULL, d1t, n_borrow = 0, pred_tps,
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
    overall  <- rwe_rmst(dta_cur, dta_ext, dta_cur_trt, n_borrow, pred_tps,
                         stderr_method)

    ## jackknife stderr
    if (stderr_method == "jk") {
        overall_theta <- overall[, 1, drop = TRUE]

        jk_theta      <- rep(0, length(overall_theta))
        for (j in seq_len(ns1)) {
            cur_jk   <- rwe_rmst(dta_cur[-j, ], dta_ext, dta_cur_trt, n_borrow,
                                 pred_tps, stderr_method)
            jk_theta <- jk_theta + (cur_jk[, 1] - overall_theta)^2
        }

        for (j in seq_len(ns1_trt)) {
            cur_jk   <- rwe_rmst(dta_cur, dta_ext, dta_cur_trt[-j, ], n_borrow,
                                 pred_tps, stderr_method)
            jk_theta <- jk_theta + (cur_jk[, 1] - overall_theta)^2
        }

        if (ns0 > 0) {
            for (j in seq_len(ns0)) {
                ext_jk   <- rwe_rmst(dta_cur, dta_ext[-j, ], dta_cur_trt,
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

#' RMST Estimation
#'
#' Estimate RMST estimates for a single PS
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
#' @return Estimation of RMST estimates at time \code{pred_tps}
#'
#'
#' @export
#'
rwe_rmst <- function(dta_cur, dta_ext, dta_cur_trt, n_borrow = 0,
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
    ## last values will be carried over for predictions
    if (is.null(pred_tps)) {
        pred_tps <- sort(unique(c(cur_surv$time[cur_surv$n.event > 0],
                                  cur_surv_trt$time[cur_surv_trt$n.event > 0])))
    }
    n_tps <- length(pred_tps)

    ## summary.survfit() only reports information for specified time points
    ## so, pred_tps should be inputted will all time points before tau.
    rst <- summary(cur_surv, time = pred_tps, extend = TRUE)
    rst_trt <- summary(cur_surv_trt, time = pred_tps, extend = TRUE)

    ## info needed for RMST
    surv_trt <- c(1, rst_trt$surv)
    surv_ctl <- c(1, rst$surv)
    time_diff <- diff(c(0, pred_tps, pred_tps[length(pred_tps)]))

    ## RMST main statistic
    area_trt <- surv_trt * time_diff
    area_ctl <- surv_ctl * time_diff
    area_d <- area_trt - area_ctl
    auc_d <- cumsum(area_d[-(n_tps + 1)])

    ## for RMST naive stderr
    if (stderr_method == "naive") {
        ## - see survival vignette page 13
        ## - \hat{\mu} = \int_0^T \hat{S}(t) dt
        ## - var(\hat{\mu}) =
        ##     \int_0^T (\int_t^T \hat{S}(u) du)^2 *
        ##              \frac{d \bar{N}(t)}{\bar{Y}(t) *
        ##                                  (\bar{Y}(t) - \bar{N}(t))}
        ## - see also similarly in survRM2:::rmst1()
        n_risk_trt <- rst_trt$n.risk
        n_risk_ctl <- rst$n.risk
        n_event_trt <- rst_trt$n.event
        n_event_ctl <- rst$n.event

        v_trt <- ifelse((n_risk_trt - n_event_trt) == 0, 0,
                        n_event_trt / (n_risk_trt * (n_risk_trt - n_event_trt)))
        v_ctl <- ifelse((n_risk_ctl - n_event_ctl) == 0, 0,
                        n_event_ctl / (n_risk_ctl * (n_risk_ctl - n_event_ctl)))

        ## This version can be very slow
        # f_getvar <- function(j, area, v) {
        #     sum(cumsum(rev(area[1:j]))^2 * rev(v[1:j]))
        # }
        # auc_var_trt <- do.call("c", lapply(1:n_tps,
        #                                    f_getvar,
        #                                    area_trt[-1], v_trt))
        # auc_var_ctl <- do.call("c", lapply(1:n_tps,
        #                                    f_getvar,
        #                                    area_ctl[-1], v_ctl))

        rev_area_trt <- rev(area_trt[-1])
        rev_area_ctl <- rev(area_ctl[-1])
        rev_v_trt <- rev(v_trt)
        rev_v_ctl <- rev(v_ctl)
        f_getvar <- function(j, rev_area, rev_v, n_tps) {
            id <- (n_tps - j + 1):n_tps
            sum(cumsum(rev_area[id])^2 * rev_v[id])
        }
        auc_var_trt <- do.call("c", lapply(1:n_tps, f_getvar,
                                           rev_area_trt, rev_v_trt, n_tps))
        auc_var_ctl <- do.call("c", lapply(1:n_tps, f_getvar,
                                           rev_area_ctl, rev_v_ctl, n_tps))

        auc_var_d <- auc_var_trt + auc_var_ctl
        auc_stderr_d <- sqrt(auc_var_d)
    } else {
        ## for none, jk, sjk, cjk, sbs, or cbs
        auc_stderr_d <- rep(NA, n_tps)
    }

    ## combine RMST estimates
    rst_rmst <- cbind(auc_d, auc_stderr_d, pred_tps)

    colnames(rst_rmst) <- c("Mean", "StdErr", "T")
    return(rst_rmst)
}

