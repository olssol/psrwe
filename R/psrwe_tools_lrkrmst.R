#' Get estimates for log-rank and RMST tests between two arms
#' for RCT augmenting control
#'
#' @noRd
#'
get_ps_lrk_rmst <- function(dta_psbor,
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
        cur_effect   <- f_stratum(cur_d1, cur_d0, cur_d1t,
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
    return(rst)
}


## Jackknife overall

#' Get JKoverall estimates for log-rank and RMST estimations between two arms
#' for RCT augmenting control
#'
#' @noRd
#'
get_ps_lrk_rmst_jkoverall <- function(dta_psbor,
                                      v_outcome     = NULL,
                                      v_event       = NULL,
                                      v_time        = NULL,
                                      f_stratum     = get_surv_stratum_lrk,
                                      f_overall_est = get_overall_est_wostderr,
                                      ...) {

    ## prepare data
    data    <- dta_psbor$data
    data    <- data[!is.na(data[["_strata_"]]), ]

    ## main estimates
    rst <- get_ps_lrk_rmst(dta_psbor, v_outcome = v_outcome,
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
        rst_jk <- get_ps_lrk_rmst(dta_psbor_jk, v_outcome = v_outcome,
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

