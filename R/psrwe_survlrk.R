#' PS-Integrated Log-Rank Test For Comparing Time-to-event Outcomes
#'
#' Log-rank tests of comparing survival outcomes up to a given time point.
#' Variance is estimated by Jack-Knife method.
#' Applies to the case when there is only one external data source and
#' two-arm RCT.
#'
#' @inheritParams psrwe_survkm
#'
#' @param v_time Column name corresponding to event time
#' @param v_event Column name corresponding to event status
#' @param pred_tp Time of interest (e.g., 1 year)
#' @param stderr_method Method for computing StdErr, see Details
#' @param ... Additional Parameters
#'
#' @details \code{stderr_method} include \code{naive} as default which
#'     mostly follows Greenwood formula (plug-in),
#'     \code{jk} using Jackknife method within each stratum, or
#'     \code{jkoverall} using Jackknife method for overall/combined estimates.
#'     Note that \code{jkoverall} may take a while longer to finish.
#'
#' @return A data frame with class name \code{PSRWE_RST_TESTANA}.
#'     It contains the test statistics of each stratum as well as the
#'     jackknife estimation. The results should be further
#'     summarized by its S3 method \code{summary}.
#'
#'
#' @examples
#' data(ex_dta_rct)
#' dta_ps_rct <- psrwe_est(ex_dta,
#'                         v_covs = paste("V", 1:7, sep = ""),
#'                         v_grp = "Group", cur_grp_level = "current",
#'                         v_arm = "Arm", ctl_arm_level = "control",
#'                         ps_method = "logistic", nstrata = 5,
#'                         stra_ctl_only = FALSE)
#' ps_bor_rct <- psrwe_borrow(dta_ps_rct, total_borrow = 30)
#' rst_lrk <- psrwe_survlrk(ps_bor_rct,
#'                          v_time = "Y_Surv",
#'                          v_event = "Status")
#' rst_lrk
#'
#' @export
#'
psrwe_survlrk <- function(dta_psbor,
                          v_time        = "time",
                          v_event       = "event",
                          pred_tp       = 1,
                          stderr_method = c("naive", "jk", "jkoverall"), 
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

    ## call estimation
    if (stderr_method %in% c("naive", "jk")) {
        rst <- get_ps_cl_lrk(dta_psbor,
                             v_event = v_event, v_time = v_time,
                             f_stratum = get_surv_stratum_lrk,
                             pred_tp = all_tps,
                             stderr_method = stderr_method,
                             ...)
    } else {
        rst <- get_ps_cl_lrk_jkoverall(dta_psbor,
                                       v_event = v_event, v_time = v_time,
                                       f_stratum = get_surv_stratum_lrk_wostderr,
                                       pred_tp = all_tps,
                                       ...)
    }

    ## return
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
get_surv_stratum_lrk <- function(d1, d0 = NULL, d1t, n_borrow = 0, pred_tp,
                                 stderr_method, ...) {

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

    ##  overall estimate
    overall  <- rwe_lrk(dta_cur, dta_ext, dta_cur_trt, n_borrow, pred_tp)

    ##jackknife stderr
    if (stderr_method == "jk") {
        overall_theta <- overall[, 1, drop = TRUE]

        jk_theta      <- rep(0, length(overall_theta))
        for (j in seq_len(ns1)) {
            cur_jk   <- rwe_lrk(dta_cur[-j, ], dta_ext, dta_cur_trt, n_borrow,
                                pred_tp)
            jk_theta <- jk_theta + (cur_jk[, 1] - overall_theta)^2
        }

        for (j in seq_len(ns1_trt)) {
            cur_jk   <- rwe_lrk(dta_cur, dta_ext, dta_cur_trt[-j, ], n_borrow,
                                pred_tp)
            jk_theta <- jk_theta + (cur_jk[, 1] - overall_theta)^2
        }

        if (ns0 > 0) {
            for (j in seq_len(ns0)) {
                ext_jk   <- rwe_lrk(dta_cur, dta_ext[-j, ], dta_cur_trt,
                                    n_borrow, pred_tp)
                jk_theta <- jk_theta + (ext_jk[, 1] - overall_theta)^2
            }
        }

        ## summary
        stderr_theta <- sqrt((ns1 + ns0 + ns1_trt - 1) / (ns1 + ns0 + ns1_trt) *
                             jk_theta)

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
#' @param dta_cur Matrix of time and event from a PS stratum in current study
#'                (control arm only)
#' @param dta_ext Matrix of time and event from a PS stratum in external data
#'                source (control arm only)
#' @param dta_cur_trt Matrix of time and event from a PS stratum in current
#'                    study (treatment arm only)
#' @param n_borrow Number of subjects to be borrowed
#' @param pred_tp Time points to be estimated
#'
#' @return Estimation of log-rank estimates at time \code{pred_tps}
#'
#'
#' @export
#'
rwe_lrk <- function(dta_cur, dta_ext, dta_cur_trt, n_borrow = 0,
                    pred_tp = 1) {

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
                        data    = cur_data,
                        weights = cur_weights)

    ## trt arm
    cur_data_trt    <- dta_cur_trt
    ns1_trt         <- nrow(dta_cur_trt)
    cur_weights_trt <- rep(1, ns1_trt)

    ## KM with stratum weights for trt arm
    colnames(cur_data_trt) <- c("time", "event")
    cur_data_trt <- data.frame(cur_data_trt)
    cur_surv_trt <- survfit(Surv(time, event) ~ 1,
                            data    = cur_data_trt,
                            weights = cur_weights_trt)

    ## summary.survfit() need to be extend to longer time points
    ## Last values will be carried over for predictions
    rst <- summary(cur_surv, time = pred_tp, extend = TRUE)
    rst_trt <- summary(cur_surv_trt, time = pred_tp, extend = TRUE)

    ## LRK
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

    ## LRK naive stderr
    std_err_d <- ifelse(n_risk <= 1, 0,
                        E_1_j * (1 - p_event) *
                        n_risk_ctl / (n_risk - 1))
    std_err_d <- sqrt(cumsum(std_err_d))

    ## Compute log-rank estimates
    rst_lrk <- cbind(mean_d, std_err_d, pred_tp)

    colnames(rst_lrk) <- c("Mean", "StdErr", "T")
    return(rst_lrk)
}


#' Get estimates for log-rank test between two arms
#'
#'
#'
#' @noRd
#'
get_ps_cl_lrk <- function(dta_psbor,
                          v_outcome     = NULL,
                          v_event       = NULL,
                          v_time        = NULL,
                          f_stratum     = get_surv_stratum_lrk,
                          f_overall_est = get_overall_est,
                          ...) {

    ## prepare data
    data    <- dta_psbor$data
    data    <- data[!is.na(data[["_strata_"]]), ]

    strata  <- levels(data[["_strata_"]])
    nstrata <- length(strata)
    borrow  <- dta_psbor$Borrow$N_Borrow

    ## estimate
    eff_theta <- NULL
    for (i in seq_len(nstrata)) {
        cur_01  <- get_cur_d(data,
                             strata[i],
                             c(v_outcome, v_time, v_event))

        cur_d1  <- cur_01$cur_d1    ## This is "cur_d1c"
        cur_d0  <- cur_01$cur_d0
        cur_d1t <- cur_01$cur_d1t

        ## effect with borrowing
        cur_effect   <- f_stratum(cur_d1t, cur_d1, cur_d0,
                                  n_borrow = borrow[i], ...)
        eff_theta <- rbind(eff_theta, cur_effect)
    }

    ## summary
    rst_effect <- f_overall_est(eff_theta, dta_psbor$Borrow$N_Current)

    ## return
    rst <-  list(Effect    = rst_effect,
                 Borrow    = dta_psbor$Borrow,
                 Total_borrow = dta_psbor$Total_borrow,
                 is_rct       = dta_psbor$is_rct)
}


## Jackknife overall

#' Get log-rank estimation for each stratum without stderr
#'
#'
#' @noRd
#'
get_surv_stratum_lrk_wostderr <- function(d1, d0 = NULL, d1t, n_borrow = 0,
                                          pred_tp, stderr_method, ...) {

    ## treatment or control only
    dta_cur <- d1
    dta_ext <- d0
    dta_cur_trt <- d1t

    ##  overall estimate
    overall  <- rwe_lrk(dta_cur, dta_ext, dta_cur_trt, n_borrow, pred_tp)
    return(overall)
}


#' Get JKoverall estimates for log-rank estimation
#'
#'
#'
#' @noRd
#'
get_ps_cl_lrk_jkoverall <- function(dta_psbor,
                                    v_outcome     = NULL,
                                    v_event       = NULL,
                                    v_time        = NULL,
                                    f_stratum     = get_surv_stratum_lrk_wostderr,
                                    f_overall_est = get_overall_est_wostderr,
                                    ...) {

    ## prepare data
    data    <- dta_psbor$data
    data    <- data[!is.na(data[["_strata_"]]), ]

    ## main estimates
    rst <- get_ps_cl_lrk(dta_psbor, v_outcome = v_outcome,
                         v_event = v_event, v_time = v_time,
                         f_stratum = f_stratum,
                         f_overall_est = f_overall_est, ...)

    ## JK overall stderr
    rstom <- rst$Effect$Overall_Estimate$Mean
    sdf <- rep(0, length(rstom))

    n_jk <- nrow(data)
    dta_psbor_jk <- dta_psbor
    for (i_jk in 1:n_jk) {
        dta_psbor_jk$data <- data[-i_jk,]
        rst_jk <- get_ps_cl_lrk(dta_psbor_jk, v_outcome = v_outcome,
                                v_event = v_event, v_time = v_time,
                                f_stratum = f_stratum,
                                f_overall_est = f_overall_est, ...)
        sdf <- sdf + (rst_jk$Effect$Overall_Estimate$Mean - rstom)^2
    }

    ## update rst
    nc_jk <- (n_jk - 1) / n_jk
    rst$Effect$Overall_Estimate$StdErr <- sqrt(sdf * nc_jk)

    ## return
    return(rst)
}

