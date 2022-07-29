#' PS-Integrated Restricted Mean Survival Time (RMST) Test For Comparing Time-to-event Outcomes
#'
#' RMST test evaluates two-arm RCT for up to a given time point.
#' Variance can be estimated by Jackknife methods.
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
#' rst_rmst <- psrwe_survrmst(ps_bor_rct,
#'                            v_time = "Y_Surv",
#'                            v_event = "Status")
#' rst_rmst
#'
#' @export
#'
psrwe_survrmst <- function(dta_psbor,
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

    ## observed
    rst_obs <- get_km_observed(data, v_time, v_event, all_tps)

    ## call estimation
    if (stderr_method %in% c("naive", "jk")) {
        rst <- get_ps_cl_rmst(dta_psbor,
                              v_event = v_event, v_time = v_time,
                              f_stratum = get_surv_stratum_rmst,
                              pred_tp = all_tps,
                              stderr_method = stderr_method,
                              ...)
    } else {
        rst <- get_ps_cl_rmst_jkoverall(dta_psbor,
                                        v_event = v_event, v_time = v_time,
                                        f_stratum = get_surv_stratum_rmst_wostderr,
                                        pred_tp = all_tps,
                                        ...)
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
get_surv_stratum_rmst <- function(d1, d0 = NULL, d1t, n_borrow = 0, pred_tp,
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
    overall  <- rwe_rmst(dta_cur, dta_ext, dta_cur_trt, n_borrow, pred_tp)

    ##jackknife stderr
    if (stderr_method == "jk") {
        overall_theta <- overall[, 1, drop = TRUE]

        jk_theta      <- rep(0, length(overall_theta))
        for (j in seq_len(ns1)) {
            cur_jk   <- rwe_rmst(dta_cur[-j, ], dta_ext, dta_cur_trt, n_borrow,
                                 pred_tp)
            jk_theta <- jk_theta + (cur_jk[, 1] - overall_theta)^2
        }

        for (j in seq_len(ns1_trt)) {
            cur_jk   <- rwe_rmst(dta_cur, dta_ext, dta_cur_trt[-j, ], n_borrow,
                                 pred_tp)
            jk_theta <- jk_theta + (cur_jk[, 1] - overall_theta)^2
        }

        if (ns0 > 0) {
            for (j in seq_len(ns0)) {
                ext_jk   <- rwe_rmst(dta_cur, dta_ext[-j, ], dta_cur_trt,
                                     n_borrow, pred_tp)
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
#' @param dta_cur Matrix of time and event from a PS stratum in current study
#'                (control arm only)
#' @param dta_ext Matrix of time and event from a PS stratum in external data
#'                source (control arm only)
#' @param dta_cur_trt Matrix of time and event from a PS stratum in current
#'                    study (treatment arm only)
#' @param n_borrow Number of subjects to be borrowed
#' @param pred_tp Time points to be estimated
#'
#' @return Estimation of RMST estimates at time \code{pred_tps}
#'
#'
#' @export
#'
rwe_rmst <- function(dta_cur, dta_ext, dta_cur_trt, n_borrow = 0,
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

    ## Info needed for RMST
    surv_trt <- rst_trt$surv
    surv_ctl <- rst$surv
    time_diff <- diff(c(0, pred_tp))
    n_risk_trt <- rst_trt$n.risk
    n_risk_ctl <- rst$n.risk
    n_event_trt <- rst_trt$n.event
    n_event_ctl <- rst$n.event

    ## RMST main statistic
    area_trt <- surv_trt * time_diff
    area_ctl <- surv_ctl * time_diff
    area_diff <- area_trt - area_ctl
    auc_diff <- cumsum(area_diff)

    ## RMST naive stderr
    ## For RMST2 naive stderr
    ## - see survival vignette page 13
    ## - \hat{\mu} = \int_0^T \hat{S}(t) dt
    ## - var(\hat{\mu}) =
    ##     \int_0^T (\int_t^T \hat{S}(u) du)^2 *
    ##              \frac{d \bar{N}(t)}{\bar{Y}(t) * (\bar{Y}(t) - \bar{N}(t))}
    ## - see also survRM2:::rmst1()
    v_trt <- ifelse((n_risk_trt - n_event_trt) == 0, 0,
                    n_event_trt / (n_risk_trt * (n_risk_trt - n_event_trt)))
    v_ctl <- ifelse((n_risk_ctl - n_event_ctl) == 0, 0,
                    n_event_ctl / (n_risk_ctl * (n_risk_ctl - n_event_ctl)))
    f_get_auc_var <- function(i.tau, area, v) {
        rev_auc <- rev(cumsum(rev(area[1:i.tau])))
        auc_var <- sum(rev_auc^2 * v[1:i.tau])
        auc_var
    }
    auc_var_trt <- do.call("c", lapply(1:length(pred_tp), f_get_auc_var,
                                       area = area_trt, v = v_trt))
    auc_var_ctl <- do.call("c", lapply(1:length(pred_tp), f_get_auc_var,
                                       area = area_ctl, v = v_ctl))
    auc_var_diff <- auc_var_trt + auc_var_ctl
    auc_stderr_diff <- sqrt(auc_var_diff)

    ## Compute RMST estimates
    rst_rmst <- cbind(auc_diff, auc_stderr_diff, pred_tp)

    colnames(rst_rmst) <- c("Mean", "StdErr", "T")
    return(rst_rmst)
}


#' Get estimates for RMST test between two arms
#'
#'
#'
#' @noRd
#'
get_ps_cl_rmst <- function(dta_psbor,
                           v_outcome     = NULL,
                           v_event       = NULL,
                           v_time        = NULL,
                           f_stratum     = get_surv_stratum_rmst,
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
    rst <-  list(Control   = NULL,
                 Treatment = NULL,
                 Effect    = rst_effect,
                 Borrow    = dta_psbor$Borrow,
                 Total_borrow = dta_psbor$Total_borrow,
                 is_rct       = dta_psbor$is_rct)
}


## Jackknife overall

#' Get RMST estimation for each stratum without stderr
#'
#'
#' @noRd
#'
get_surv_stratum_rmst_wostderr <- function(d1, d0 = NULL, d1t, n_borrow = 0,
                                           pred_tp, stderr_method, ...) {

    ## treatment or control only
    dta_cur <- d1
    dta_ext <- d0
    dta_cur_trt <- d1t

    ##  overall estimate
    overall  <- rwe_rmst(dta_cur, dta_ext, dta_cur_trt, n_borrow, pred_tp)
    return(overall)
}


#' Get JKoverall estimates for RMST estimation
#'
#'
#'
#' @noRd
#'
get_ps_cl_rmst_jkoverall <- function(dta_psbor,
                                     v_outcome     = NULL,
                                     v_event       = NULL,
                                     v_time        = NULL,
                                     f_stratum     = get_surv_stratum_rmst_wostderr,
                                     f_overall_est = get_overall_est_wostderr,
                                     ...) {

    ## prepare data
    data    <- dta_psbor$data
    data    <- data[!is.na(data[["_strata_"]]), ]

    ## main estimates
    rst <- get_ps_cl_rmst(dta_psbor, v_outcome = v_outcome,
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
        rst_jk <- get_ps_cl_rmst(dta_psbor_jk, v_outcome = v_outcome,
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

