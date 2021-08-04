#' PS-Integrated Kaplan-Meier Estimation Plot
#'
#' Plot the mean of a survival outcome at all distinctive time points based on
#' PS-integrated Kaplan-Meier approach. Variance is estimated by usual Greenwood
#' method. Applies to the case when there is only one external data source.
#'
#' @inheritParams rwe_ps_powerp
#'
#' @param v_time Column name corresponding to event time
#' @param v_event Column name corresponding to event status
#'
#' @param ... Additional Parameters.
#'
#' @return A data frame with class name \code{RWE_PS_RST}. It contains the
#'     composite estimation of the mean for overall and each stratum.
#'
#'
#' @examples
#' \donttest{
#' data(ex_dta)
#' dta_ps <- rwe_ps_est(ex_dta,
#'        v_covs = paste("V", 1:7, sep = ""),
#'        v_grp = "Group",
#'        cur_grp_level = "current")
#' ps_borrow <- rwe_ps_borrow(total_borrow = 30, dta_ps)
#' rst       <- rwe_ps_survkmplot(ps_borrow, v_time = "Y_time",
#'                                v_event = "Status")}
#'
#' @export
#'
rwe_ps_survkmplot <- function(dta_psbor,
                              v_time     = "time",
                              v_event    = "event",
                              ...) {

    ## check
    stopifnot(inherits(dta_psbor,
                       what = get_rwe_class("PSDIST")))

    stopifnot(all(c(v_event, v_time) %in%
                  colnames(dta_psbor$data)))

    pred_tp <- sort(unique(v_time))

    ## observed
    rst_obs <- get_km_observed(dta_psbor$data, v_time, v_event, pred_tp)

    ## call estimation
    rst <- get_ps_cl_kmplot(dta_psbor, v_event = v_event, v_time = v_time,
                            f_stratum = get_surv_stratum, pred_tp = pred_tp,
                            ...)

    ## return
    rst$Observed <- rst_obs
    rst$pred_tp  <- pred_tp
    rst$Method   <- "ps_km"
    class(rst)   <- get_rwe_class("ANARSTPLT")
    rst
}


#' Get estimates for survival at all distinctive time points
#'
#'
#'
#' @noRd
#'
get_ps_cl_kmplot <- function(dta_psbor,
                             v_outcome  = NULL,
                             v_event    = NULL,
                             v_time     = NULL,
                             f_stratum  = get_cl_stratum,
                             ...) {
### TBD

    ## prepare data
    is_rct  <- dta_psbor$is_rct
    data    <- dta_psbor$data
    data    <- data[!is.na(data[["_strata_"]]), ]

    strata  <- levels(data[["_strata_"]])
    nstrata <- length(strata)
    borrow  <- dta_psbor$Borrow$N_Borrow

    ## estimate
    ctl_theta <- NULL
    trt_theta <- NULL

    for (i in seq_len(nstrata)) {
        cur_01  <- get_cur_d(data,
                             strata[i],
                             c(v_outcome, v_time, v_event))

        cur_d1  <- cur_01$cur_d1
        cur_d0  <- cur_01$cur_d0
        cur_d1t <- cur_01$cur_d1t

        ## control with borrowing
        cur_ctl   <- f_stratum(cur_d1, cur_d0, n_borrow = borrow[i], ...)
        ctl_theta <- rbind(ctl_theta, cur_ctl)
        if (is_rct) {
            cur_trt   <- f_stratum(cur_d1t, ...)
            trt_theta <- rbind(trt_theta, cur_trt)
        }
    }

    ## summary
    rst_trt    <- NULL
    rst_effect <- NULL
    if (is_rct) {
        rst_trt    <- get_overall_est(trt_theta[, 1],
                                      trt_theta[, 2],
                                      dta_psbor$Borrow$N_Cur_TRT)
        rst_effect <- get_overall_est(trt_theta[, 1] - ctl_theta[, 1],
                                      sqrt(trt_theta[, 2] + ctl_theta[, 2]),
                                      dta_psbor$Borrow$N_Current)
        n_ctl      <- dta_psbor$Borrow$N_Cur_CTL
    } else {
        n_ctl      <- dta_psbor$Borrow$N_Current
    }
    rst_ctl <- get_overall_est(ctl_theta[, 1], ctl_theta[, 2], n_ctl)

    ## return
    rst <-  list(Control   = rst_ctl,
                 Treatment = rst_trt,
                 Effect    = rst_effect,
                 Borrow    = dta_psbor$Borrow,
                 Total_borrow = dta_psbor$Total_borrow,
                 is_rct       = is_rct)
}

