#' PS-Integrated Kaplan-Meier Estimation
#'
#' Estimate the mean of a survival outcome at a given time point based on
#' PS-integrated Kaplan-Meier approach. Variance can be estimated by the
#' Jackknife methods. Apply to the case when there is only one external
#' data source.
#'
#' @inheritParams psrwe_powerp
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
#' @return A data frame with class name \code{PSRWE_RST}. It contains the
#'     composite estimation of the mean for each stratum as well as the
#'     Jackknife estimation. The results can be further
#'     summarized by its S3 method \code{summary}.
#'     The results can also be analyzed by \code{psrwe_outana} for outcome
#'     analysis and inference.
#'
#'
#' @examples
#' data(ex_dta)
#' dta_ps <- psrwe_est(ex_dta,
#'        v_covs = paste("V", 1:7, sep = ""),
#'        v_grp = "Group",
#'        cur_grp_level = "current")
#' ps_borrow <- psrwe_borrow(total_borrow = 30, dta_ps)
#' rst       <- psrwe_survkm(ps_borrow,
#'                           pred_tp = 365,
#'                           v_time = "Y_Surv",
#'                           v_event = "Status")
#' rst
#'
#' @export
#'
psrwe_survkm <- function(dta_psbor, pred_tp,
                         v_time     = "time",
                         v_event    = "event",
                         stderr_method = c("naive", "jk", "sjk", "cjk",
                                           "sbs", "cbs", "none"), 
                         n_bootstrap = 200,
                         ...) {

    ## check
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
    if (stderr_method[1] %in% c("naive", "jk", "none")) {
        rst <- get_ps_cl_km(dta_psbor,
                            v_event = v_event, v_time = v_time,
                            f_stratum = get_surv_stratum,
                            pred_tp = all_tps,
                            stderr_method = stderr_method,
                            ...)
    } else if (stderr_method[1] %in% c("sjk")) {
        rst <- get_ps_cl_km_sjk(dta_psbor,
                                v_event = v_event, v_time = v_time,
                                f_stratum = get_surv_stratum,
                                pred_tp = all_tps,
                                stderr_method = "none",
                                ...)
    } else if (stderr_method[1] %in% c("cjk")) {
        rst <- get_ps_cl_km_cjk(dta_psbor,
                                v_event = v_event, v_time = v_time,
                                f_stratum = get_surv_stratum,
                                pred_tp = all_tps,
                                stderr_method = "none",
                                ...)
    } else if (stderr_method[1] %in% c("sbs")) {
        rst <- get_ps_cl_km_sbs(dta_psbor,
                                v_event = v_event, v_time = v_time,
                                f_stratum = get_surv_stratum,
                                pred_tp = all_tps,
                                stderr_method = "none",
                                n_bootstrap = n_bootstrap,
                                ...)
    } else if (stderr_method[1] %in% c("cbs")) {
        rst <- get_ps_cl_km_cbs(dta_psbor,
                                v_event = v_event, v_time = v_time,
                                f_stratum = get_surv_stratum,
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
    rst$Method   <- "ps_km"
    rst$Outcome_type <- "tte"
    class(rst)   <- get_rwe_class("ANARST")
    return(rst)
}

#' Get surv estimation for each stratum
#'
#'
#' @noRd
#'
get_surv_stratum <- function(d1, d0 = NULL, n_borrow = 0, pred_tps,
                             stderr_method = "jk", ...) {

    ## treatment or control only
    dta_cur <- d1
    dta_ext <- d0

    ## overall estimate
    overall  <- rwe_km(dta_cur,
                       dta_ext       = dta_ext,
                       n_borrow      = n_borrow,
                       pred_tps      = pred_tps,
                       stderr_method = stderr_method)

    ## jackknife stderr
    if (stderr_method == "jk") {
        ns1     <- nrow(dta_cur)
        if (is.null(d0)) {
            ns0 <- 0
        } else {
            ns0 <- nrow(dta_ext)
        }

        overall_theta <- overall[, 1, drop = TRUE]

        jk_theta      <- rep(0, length(overall_theta))
        for (j in seq_len(ns1)) {
            cur_jk   <- rwe_km(dta_cur[-j, ],
                               dta_ext       = dta_ext,
                               n_borrow      = n_borrow,
                               pred_tps      = pred_tps,
                               stderr_method = stderr_method)
            jk_theta <- jk_theta + (cur_jk[, 1] - overall_theta)^2
        }

        if (ns0 > 0) {
            for (j in seq_len(ns0)) {
                ext_jk   <- rwe_km(dta_cur,
                                   dta_ext       = dta_ext[-j, ],
                                   n_borrow      = n_borrow,
                                   pred_tps      = pred_tps,
                                   stderr_method = stderr_method)
                jk_theta <- jk_theta + (ext_jk[, 1] - overall_theta)^2
            }
        }

        ## summary
        stderr_theta <- sqrt((ns1 + ns0 - 1) / (ns1 + ns0) * jk_theta)

        overall[, 2] <- stderr_theta
    }

    return(overall)
}

#' Kaplan-Meier Estimation
#'
#' Estimate survival probability based on Kaplan-Meier estimator for a single PS
#' stratum
#'
#'
#' @param dta_cur Matrix of time and event from a PS stratum in the current
#'     study
#' @param dta_ext Matrix of time and event from a PS stratum in the external
#'     data source
#' @param n_borrow Number of subjects to be borrowed
#' @param pred_tps Time points to be estimated (unique and sorted)
#' @param stderr_method Method for computing StdErr (available for naive only)
#'
#' @return Estimation of survival probabilities at time \code{pred_tps}
#'
#'
#' @export
#'
rwe_km <- function(dta_cur, dta_ext = NULL, n_borrow = 0, pred_tps = NULL,
                   stderr_method = "naive") {

    ## current control and external control if available
    cur_data    <- dta_cur
    ns1         <- nrow(dta_cur)
    cur_weights <- rep(1, ns1)

    if (!is.null(dta_ext) & n_borrow > 0) {
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

    ## summary.survfit() need to be extend to longer time points
    ## Last values will be carried over for predictions
    if (is.null(pred_tps)) {
        pred_tps <- sort(unique(c(cur_surv$time[cur_surv$n.event > 0])))
    }

    ## summary.survfit() need to be extend to longer time points
    ## Last values will be carried over for predictions
    rst <- summary(cur_surv, time = pred_tps, extend = TRUE)

    ## For KM naive stderr
    if (stderr_method == "naive") {
        rst <- cbind(rst$surv, rst$std.err, pred_tps)
    } else {
        ## For none, jk, sjk, cjk, sbs, or cbs
        rst <- cbind(rst$surv, rep(NA, length(rst$std.err)), pred_tps)
    }

    colnames(rst) <- c("Mean", "StdErr", "T")
    return(rst)
}

#' Get observed KM
#'
#' @noRd
#'
get_km_observed <- function(dta, v_time, v_event, pred_tp) {
    rst <- NULL
    rst_overall <- NULL
    for (g in unique(dta[["_grp_"]])) {
        for (a in unique(dta[["_arm_"]])) {
            cur_d <- dta %>%
                dplyr::filter(g == `_grp_` &
                              a == `_arm_`)

            if (0 == nrow(cur_d))
                next

            for (s in levels(dta[["_strata_"]])) {
                cur_s <- cur_d %>%
                    dplyr::filter(s == `_strata_`)

                if (0 == nrow(cur_s))
                    next

                est <- rwe_km(cur_s[, c(v_time, v_event), drop = F],
                              pred_tps = pred_tp)

                rst <- rbind(rst,
                             data.frame(Group   = g,
                                        Arm     = a,
                                        Stratum = s,
                                        N       = nrow(cur_s),
                                        Mean    = est[, 1],
                                        StdErr  = est[, 2],
                                        T       = est[, 3]))
            }

            est <- rwe_km(cur_d[, c(v_time, v_event)],
                          pred_tps = pred_tp)
            rst_overall <- rbind(rst_overall,
                                 data.frame(Group   = g,
                                            Arm     = a,
                                            Stratum = "Overall",
                                            N       = nrow(cur_d),
                                            Mean    = est[, 1],
                                            StdErr  = est[, 2],
                                            T       = est[, 3]))
        }
    }

    rst <- rbind(rst, rst_overall)
    return(rst)
}
