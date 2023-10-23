#' PS-Integrated Composite Likelihood Estimation (WATT)
#'
#' Estimate the mean of the outcome based on PS-integrated composite likelihood
#' approach with weights of ATT (WATT).
#' Variance is estimated by Jack-Knife method. Applies to the case
#' when there is only one external data source.
#'
#' @inheritParams psrwe_powerp
#'
#' @param stderr_method Method for computing StdErr, see Details
#' @param n_bootstrap Number of bootstrap samples (for bootstrap stderr)
#' @param ... Parameters for \code{rwe_cl_watt}
#'
#' @details \code{stderr_method} include \code{jk} as default
#'     using Jackknife method within each stratum,
#'     \code{sjk} for simple Jackknife method for combined estimates
#'     such as point estimates in single arm or treatment effects in RCT, or
#'     \code{cjk} for complex Jackknife method including refitting PS model,
#'     matching, trimming, calculating borrowing parameters, and
#'     combining overall estimates.
#'     Note that \code{sjk} may take a while longer to finish and
#'     \code{cjk} will take even much longer to finish.
#'     The \code{sbs} and \code{cbs} is for simple and complex Bootstrap
#'     methods.
#'
#' @return A data frame with class name \code{PSRWE_RST}. It contains the
#'     composite estimation of the mean for each stratum as well as the
#'     jackknife estimation for each subject. The results can be further
#'     summarized by its S3 method \code{summary}.
#'     The results can be also analyzed by \code{psrwe_outana} for outcome
#'     analysis and inference.
#'
#' @examples
#' data(ex_dta)
#' dta_ps <- psrwe_est(ex_dta,
#'        v_covs = paste("V", 1:7, sep = ""),
#'        v_grp = "Group",
#'        cur_grp_level = "current",
#'        nstrata = 1)
#' ps_borrow <- psrwe_borrow(total_borrow = 30, dta_ps)
#' rst       <- psrwe_compl_watt(ps_borrow, v_outcome = "Y_Bin")
#' rst
#'
#' @export
#'
psrwe_compl_watt <- function(dta_psbor, v_outcome = "Y",
                             outcome_type = c("continuous", "binary"),
                             stderr_method = c("jk", "sjk", "cjk",
                                               "sbs", "cbs", "none"), 
                             n_bootstrap = 200,
                             ...) {

    ## check
    stopifnot(inherits(dta_psbor,
                       what = get_rwe_class("PSDIST")))

    outcome_type <- match.arg(outcome_type)

    stopifnot(v_outcome %in% colnames(dta_psbor$data))

    stderr_method <- match.arg(stderr_method)

    ## observed
    rst_obs <- get_observed(dta_psbor$data, v_outcome)

    ## call estimation
    f_get_ps_cl_km_watt <- switch(stderr_method[1],
                                  jk = get_ps_cl_km_watt,
                                  sjk = get_ps_cl_km_watt_sjk,
                                  cjk = get_ps_cl_km_watt_cjk,
                                  sbs = get_ps_cl_km_watt_sbs,
                                  cbs = get_ps_cl_km_watt_cbs,
                                  none = get_ps_cl_km_watt_none)
    rst <- f_get_ps_cl_km_watt(dta_psbor, v_outcome = v_outcome,
                               outcome_type = outcome_type,
                               f_stratum = get_cl_stratum_watt,
                               ...)

    ## return
    rst$Observed      <- rst_obs
    rst$stderr_method <- stderr_method
    rst$Method        <- "ps_cl"
    rst$Method_weight <- "WATT"
    rst$Outcome_type  <- outcome_type
    class(rst)   <- get_rwe_class("ANARST")

    rst
}


#' Composite Likelihood Estimation (WATT)
#'
#' Estimate parameter of interest based composite likelihood for a single PS
#' stratum with weights of ATT (WATT).
#'
#' @inheritParams psrwe_powerp
#'
#' @param dta_cur Vector of outcome from a PS stratum in current study
#' @param dta_ext Vector of outcome from a PS stratum in external data source
#' @param n_borrow Number of subjects to be borrowed
#' @param dta_ext_watt_di Weights of ATT for subjects in external data source
#' @param equal_sd Boolean. whether sd is the same between the current study and
#'     external data source
#'
#' @return Maximum composite likelihood estimator of the mean
#'
#' @examples
#' x <- rnorm(100,  mean = 0, sd = 1)
#' y <- rnorm(1000, mean = 1, sd = 2)
#' rwe_cl_watt(x, y, n_borrow = 20, equal_sd = FALSE)
#'
#' @export
#'
rwe_cl_watt <- function(dta_cur, dta_ext, n_borrow = 0,
                        dta_ext_watt_di = NULL,
                        outcome_type = c("continuous", "binary"),
                        equal_sd = TRUE) {

    f_ll <- function(pars) {
        theta  <- pars[1]
        sig2_1 <- pars[2]
        sig2_0 <- pars[3]

        ll <- - n1 * log(sig2_1) / 2
        ll <- ll - n1 * mean((dta_cur - theta)^2) / 2 / sig2_1
        ll <- ll - n_borrow * log(sig2_0) / 2
        # ll <- ll - n_borrow * mean((dta_ext - theta)^2) / 2 / sig2_0
        ll <- ll - n_borrow * sum((dta_ext - theta)^2 * dta_ex_watt_di) / 2 / sig2_0

        ll
    }

    f_gradient <- function(pars) {
        theta  <- pars[1]
        sig2_1 <- pars[2]
        sig2_0 <- pars[3]

        g <- numeric(length(pars))

        ## d logl / d theta
        g[1] <- n1 / sig2_1 * (mean(dta_cur) - theta) +
            n_borrow / sig2_0 * (mean(dta_ext) - theta)

        ## d logl / d sig2.1
        g[2] <- - n1 / 2 / sig2_1 +
            n1 * mean((dta_cur - theta)^2) / 2 / sig2_1 / sig2_1
        ## d logl / d sig2.0
        g[3] <- - n_borrow / 2 / sig2_0 +
            n_borrow * mean((dta_cur - theta)^2) / 2 / sig2_0 / sig2_0

        return(g)
    }


    type <- match.arg(outcome_type)
    n1   <- length(dta_cur)

    ## ignore external data
    if (0 == n_borrow) {
        ## placeholder
        dta_ext  <- dta_cur
        equal_sd <- TRUE
    }

    # init_theta <- (n1 / (n1 + n_borrow)) * mean(dta_cur) +
    #     (n_borrow / (n1 + n_borrow)) * mean(dta_ext)
    if (is.null(dta_ext_watt_di)) {
        dta_ext_watt_di <- rep(1 / length(dta_ext), length(dta_ext))
    }
    init_theta <- (n1 / (n1 + n_borrow)) * mean(dta_cur) +
                  (n_borrow / (n1 + n_borrow)) * sum(dta_ext * dta_ext_watt_di)

    if (("continuous" == type & equal_sd) |
        "binary" == type) {
        rst <- init_theta
    } else {
        init_sig2_1 <- mean((dta_cur - init_theta)^2)
        init_sig2_0 <- mean((dta_ext - init_theta)^2)
        rst         <- optim(c(init_theta, init_sig2_1, init_sig2_0),
                             method  = "L-BFGS-B",
                             fn      = f_ll,
                             lower   = c(-Inf, 1e-6, 1e-6),
                             upper   = rep(Inf, 3),
                             control = list(fnscale = -1))$par[1];
    }

    rst
}


#' Get CL estimation for each stratum (WATT)
#'
#'
#' @noRd
#'
get_cl_stratum_watt <- function(d1, d0 = NULL, n_borrow = 0, outcome_type,
                                stderr_method = "jk",
                                d0_watt_di = NULL,
                                ...) {

    ## treatment or control only
    dta_cur <- d1
    ns1     <- length(dta_cur)
    if (0 == n_borrow | is.null(d0)) {
        theta    <- mean(dta_cur)
        stderr_theta <- switch(outcome_type,
                               continuous = sd(dta_cur) / sqrt(ns1),
                               binary     = sqrt(theta * (1 - theta) / ns1))

        return(c(theta, stderr_theta))
    }

    ## overall ps-cl
    dta_ext <- d0
    ns0     <- length(dta_ext)
    dta_ext_watt_di <- d0_watt_di

    ## overall estimate
    overall_theta  <- rwe_cl_watt(dta_cur, dta_ext, n_borrow = n_borrow,
                                  dta_ext_watt_di = dta_ext_watt_di,
                                  ...)

    ## jackknife stderr
    if (stderr_method == "jk") {
        jk_theta <- NULL
        for (j in seq_len(ns1)) {
            cur_jk   <- rwe_cl_watt(dta_cur[-j], dta_ext, n_borrow = n_borrow,
                                    dta_ext_watt_di = dta_ext_watt_di,
                                    ...)
            jk_theta <- c(jk_theta, cur_jk)
        }

        if (ns0 > 0) {
            for (j in seq_len(ns0)) {
                cur_jk <- rwe_cl_watt(dta_cur, dta_ext[-j], n_borrow = n_borrow,
                                      dta_ext_watt_di = dta_ext_watt_di[-j],
                                      ...)
                jk_theta <- c(jk_theta, cur_jk)
            }
        }

        stderr_theta <- get_jk_sd(overall_theta, jk_theta)
    } else {
        stderr_theta <- NA
    }

    ## summary
    return(c(overall_theta, stderr_theta))
}


#' Get estimates for composite likelihood and survival (WATT)
#'
#'
#'
#' @noRd
#'
get_ps_cl_km_watt <- function(dta_psbor,
                              v_outcome     = NULL,
                              v_event       = NULL,
                              v_time        = NULL,
                              f_stratum     = get_cl_stratum_watt,
                              f_overall_est = get_overall_est,
                              ...) {

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

        cur_01_e  <- get_cur_d(data, strata[i], "_ps_")
        cur_d0_e  <- cur_01_e$cur_d0
        cur_d0_watt  <- cur_d0_e / (1 - cur_d0_e)
        cur_d0_watt_di  <- cur_d0_watt / sum(cur_d0_watt)

        ## control with borrowing
        cur_ctl   <- f_stratum(cur_d1, cur_d0, n_borrow = borrow[i],
                               d0_watt_di = cur_d0_watt_di,
                               ...)
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
        rst_trt    <- f_overall_est(trt_theta, dta_psbor$Borrow$N_Cur_TRT)
        rst_effect <- f_overall_est(trt_theta, dta_psbor$Borrow$N_Current,
                                    ctl_theta)
        n_ctl      <- dta_psbor$Borrow$N_Cur_CTL
    } else {
        n_ctl      <- dta_psbor$Borrow$N_Current
    }
    rst_ctl <- f_overall_est(ctl_theta, n_ctl)

    ## return
    rst <-  list(Control   = rst_ctl,
                 Treatment = rst_trt,
                 Effect    = rst_effect,
                 Borrow    = dta_psbor$Borrow,
                 Total_borrow = dta_psbor$Total_borrow,
                 is_rct       = is_rct)
    return(rst)
}


#' Get estimates for composite likelihood and survival (WATT) skip stderr
#'
#'
#'
#' @noRd
#'
get_ps_cl_km_watt_none <- function(dta_psbor,
                                   v_outcome     = NULL,
                                   v_event       = NULL,
                                   v_time        = NULL,
                                   f_stratum     = get_cl_stratum_watt,
                                   f_overall_est = get_overall_est,
                                   ...) {
    get_ps_cl_km_watt(dta_psbor,
                      v_outcome     = v_outcome,
                      v_event       = v_event,
                      v_timet       = v_time,
                      v_stratum     = v_stragum,
                      f_overall_est = f_overall_est,
                      stderr_method = "none",
                      ...)
}
