#' Get JKoverall estimates for composite likelihood and survival
#'
#' @noRd
#'
get_ps_cl_km_jkoverall <- function(dta_psbor,
                                   v_outcome     = NULL,
                                   v_event       = NULL,
                                   v_time        = NULL,
                                   f_stratum     = get_cl_stratum,
                                   f_overall_est = get_overall_est,
                                   f_stratum_wostderr     = get_cl_stratum_wostderr,
                                   f_overall_est_wostderr = get_overall_est_wostderr,
                                   ...) {

    ## prepare data
    is_rct  <- dta_psbor$is_rct
    data    <- dta_psbor$data
    data    <- data[!is.na(data[["_strata_"]]), ]

    ## main estimates
    rst <- get_ps_cl_km(dta_psbor, v_outcome = v_outcome,
                        v_event = v_event, v_time = v_time,
                        f_stratum = f_stratum,
                        f_overall_est = f_overall_est, ...)

    ## JK overall stderr
    rstom_ctl <- rst$Control$Overall_Estimate$Mean
    sdf_ctl <- rep(0, length(rstom_ctl))

    if (is_rct) {
        rstom_trt <- rst$Treatment$Overall_Estimate$Mean
        rstom_eff <- rst$Effect$Overall_Estimate$Mean
        sdf_trt <- rep(0, length(rstom_trt))
        sdf_eff <- rep(0, length(rstom_eff))
    }

    n_jk <- nrow(data)
    dta_psbor_jk <- dta_psbor
    for (i_jk in 1:n_jk) {
        dta_psbor_jk$data <- data[-i_jk,]
        rst_jk <- get_ps_cl_km(dta_psbor_jk, v_outcome = v_outcome,
                               v_event = v_event, v_time = v_time,
                               f_stratum = f_stratum_wostderr,
                               f_overall_est = f_overall_est_wostderr, ...)

        sdf_ctl <- sdf_ctl + (rst_jk$Control$Overall_Estimate$Mean -
                              rstom_ctl)^2

        if (is_rct) {
            sdf_trt <- sdf_trt + (rst_jk$Treatment$Overall_Estimate$Mean -
                                  rstom_trt)^2
            sdf_eff <- sdf_eff + (rst_jk$Effect$Overall_Estimate$Mean -
                                  rstom_eff)^2
        }
    }

    ## update rst
    nc_jk <- (n_jk - 1) / n_jk
    rst$Control$Overall_Estimate$StdErr <- sqrt(sdf_ctl * nc_jk)

    if (is_rct) {
        rst$Treatment$Overall_Estimate$StdErr <- sqrt(sdf_trt * nc_jk)
        rst$Effect$Overall_Estimate$StdErr <- sqrt(sdf_eff * nc_jk)
    }

    ## return
    return(rst)
}


#' Summarize overall theta without stderr
#'
#'
#' @noRd
#'
get_overall_est_wostderr <- function(ts1, weights, ts2 = NULL) {

    if (is.null(ts2)) {
        theta0 <- ts1[, 1]
    } else {
        theta0 <- ts1[, 1] - ts2[, 1]
    }

    ws         <- weights / sum(weights)
    nstrata    <- length(ws)
    theta      <- matrix(theta0, ncol = nstrata)

    overall    <- as.vector(theta %*% ws)

    ## stratum est
    s_est <- data.frame(Mean   = theta0,
                        StdErr = NA)
    o_est <- data.frame(Mean   = overall,
                        StdErr = NA)

    if (ncol(ts1) > 2) {
        prept <- matrix(ts1[, 3], ncol = nstrata)
        s_est <- cbind(s_est, T = prept[, 1],
                       Stratum = rep(1:nstrata, each = nrow(theta)))
        o_est <- cbind(o_est, T = prept[, 1])
    }

    list(Stratum_Estimate = s_est,
         Overall_Estimate = o_est)
}


#' Get CL estimation for each stratum without stderr
#'
#'
#' @noRd
#'
get_cl_stratum_wostderr <- function(d1, d0 = NULL, n_borrow = 0, outcome_type, ...) {

    ## treatment or control only
    dta_cur <- d1
    ns1     <- length(dta_cur)
    if (0 == n_borrow | is.null(d0)) {
        theta    <- mean(dta_cur)
        return(c(theta, NA))
    }

    ## overall ps-cl
    dta_ext <- d0
    ns0     <- length(dta_ext)

    ##  overall estimate
    overall_theta  <- rwe_cl(dta_cur, dta_ext, n_borrow, ...)

    ##jackknife
    jk_theta <- NULL
    for (j in seq_len(ns1)) {
        cur_jk   <- rwe_cl(dta_cur[-j], dta_ext, n_borrow, ...)
        jk_theta <- c(jk_theta, cur_jk)
    }

    if (ns0 > 0) {
        for (j in seq_len(ns0)) {
            cur_jk <- rwe_cl(dta_cur, dta_ext[-j], n_borrow, ...)
            jk_theta <- c(jk_theta, cur_jk)
        }
    }

    ## summary
    return(c(overall_theta, NA))
}


#' Get surv estimation for each stratum without stderr
#'
#'
#' @noRd
#'
get_surv_stratum_wostderr <- function(d1, d0 = NULL, n_borrow = 0, pred_tp,
                                      stderr_method, ...) {

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
}


