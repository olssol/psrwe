#' PS-Integrated Kaplan-Meier Estimation
#'
#' Estimate the mean of a survival outcome at a given time point based on
#' PS-integrated Kaplan-Meier approach. Variance is estimated by Jack-Knife
#' method. Applies to the case when there is only one external data source.
#'
#' @inheritParams rwe_ps_powerp
#'
#' @param v_event Column name corresponding to event time
#' @param v_censored Column name corresponding to censored status
#' @param pred_tp Time of interest (e.g., 1 year)
#'
#' @param ... Additional Parameters.
#'
#' @return A data frame with class name \code{RWE_PS_RST}. It contains the
#'     composite estimation of the mean for each stratum as well as the
#'     jackknife estimation for each subject. The results should be further
#'     summarized by its S3 method \code{summary}.
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
#' rst       <- rwe_ps_cl(ps_borrow, v_outcome = "Y_Con")}
#'
#' @export
#'
rwe_ps_survkm <- function(dta_psbor,
                          v_time     = "time",
                          v_event    = "event",
                          pred_tp    = 1,
                          ...) {

    ## check
    stopifnot(inherits(dta_psbor,
                       what = get_rwe_class("PSDIST")))

    stopifnot(all(c(v_event, v_time) %in%
                  colnames(dta_psbor$data)))

    stopifnot(is.numeric(pred_tp) & 1 == length(pred_tp))

    ## observed
    rst_obs <- get_km_observed(dta_psbor$data, v_time, v_event, pred_tp)

    ## call estimation
    rst <- get_ps_cl_km(dta_psbor, v_event = v_event, v_time = v_time,
                        f_stratum = get_surv_stratum, pred_tp = pred_tp,
                        ...)

    ## return
    rst$Observed <- rst_obs
    rst$pred_tp  <- pred_tp
    rst$Method   <- "ps_km"
    class(rst)   <- get_rwe_class("ANARST")
    rst
}

#' Get estimation for each stratum
#'
#'
#' @noRd
#'
get_surv_stratum <- function(d1, d0 = NULL, n_borrow = 0, pred_tp, ...) {

    ## treatment or control only
    dta_cur <- d1
    dta_ext <- d0
    ns1     <- nrow(dta_cur)

    if (is.null(d0)) {
        ns0 <- 0
    } else {
        ns0 <- nrow(dta_ext)
    }

    ##  overall estimate
    overall  <- rwe_km(dta_cur, dta_ext, n_borrow, pred_tp)

    return(overall)

    if (0) {
        ## instead of using jackknife, it is using the std.err from the surv
        ## function with the same weights
        if (0 == ns0) {
            return(overall)
        }

        ##jackknife
        overall_theta <- overall[, 1, drop = TRUE]
        overall_sd    <- overall[, 2]

        jk_theta      <- NULL
        for (j in seq_len(ns1)) {
            cur_jk   <- rwe_km(dta_cur[-j, ], dta_ext, n_borrow, pred_tp)
            jk_theta <- rbind(jk_theta, cur_jk[, 1])
        }

        if (ns0 > 0) {
            for (j in seq_len(ns0)) {
                cur_jk   <- rwe_km(dta_cur, dta_ext[-j, ], n_borrow, pred_tp)
                jk_theta <- rbind(jk_theta, cur_jk[, 1])
            }
        }

        ## summary
        sd_theta <- apply(rbind(overall_theta, jk_theta),
                          2,
                          function(x) get_jk_sd(x[1], x[-1]))

        return(cbind(overall_theta, sd_theta))
    }
}

#' Kaplan-Meier Estimation
#'
#' Estimate survival probability based on Kaplan-Meier estimator for a single PS
#' stratum
#'
#'
#' @param dta_cur Matrix of time and event from a PS stratum in current study
#' @param dta_ext Matrix of time and event from a PS stratum in external data source
#' @param n_borrow Number of subjects to be borrowed
#' @param pred_tps Time points to be estimated
#'
#' @return Estimation of survival probabilities at time \code{pred_tps}
#'
#'
#' @export
#'
rwe_km <- function(dta_cur, dta_ext = NULL, n_borrow = 0, pred_tp = 1) {

    cur_data    <- dta_cur
    ns1         <- nrow(dta_cur)
    cur_weights <- rep(1, ns1)

    if (!is.null(dta_ext) & n_borrow > 0) {
        ns0         <- nrow(dta_ext)
        cur_data    <- rbind(cur_data, dta_ext)
        cur_weights <- c(cur_weights,
                         rep(n_borrow / ns0, ns0))
    }

    colnames(cur_data) <- c("time", "event")
    cur_data  <- data.frame(cur_data)
    cur_surv <- survfit(Surv(time, event) ~ 1,
                        data    = cur_data,
                        weights = cur_weights)

    rst <- summary(cur_surv, time = pred_tp)
    rst <- cbind(rst$surv, rst$std.err)

    if (nrow(rst) < length(pred_tp)) {
        inx <- 1:nrow(rst)
        inx <- c(inx,
                 rep(nrow(rst), length(pred_tp) - nrow(rst)))

        rst <- rst[inx, ]
    }

    rst
}

#' Get observed KM
#'
#' @noRd
#'
get_km_observed <- function(dta, v_time, v_event, pred_tp) {

    rst <- NULL
    for (g in unique(dta[["_grp_"]])) {
        for (a in unique(dta[["_arm_"]])) {
            cur_d <- dta %>%
                dplyr::filter(g == `_grp_` &
                              a == `_arm_`)

            if (0 == nrow(cur_d))
                next

            est <- rwe_km(cur_d[, c(v_time, v_event)], pred_tp = pred_tp)
            rst <- rbind(rst,
                         data.frame(Group   = g,
                                    Arm     = a,
                                    Stratum = "Overall",
                                    N       = nrow(cur_d),
                                    Mean    = est[1],
                                    SD      = est[2]))

            for (s in levels(dta[["_strata_"]])) {
                cur_s <- cur_d %>%
                    dplyr::filter(s == `_strata_`)

                if (0 == nrow(cur_s))
                    next

                est <- rwe_km(cur_s[, c(v_time, v_event), drop = F],
                              pred_tp = pred_tp)

                rst <- rbind(rst,
                             data.frame(Group   = g,
                                        Arm     = a,
                                        Stratum = s,
                                        N       = nrow(cur_s),
                                        Mean    = est[1],
                                        SD      = est[2]))
            }
        }
    }

    rst

}
