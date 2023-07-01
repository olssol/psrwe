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


## Simple Jackknife

#' Get JKoverall estimates for log-rank and RMST estimations between two arms
#' for RCT augmenting control
#'
#' @noRd
#'
get_ps_lrk_rmst_sjk <- function(dta_psbor,
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
    rst_om <- rst$Effect$Overall_Estimate$Mean
    sdf_om <- rep(0, length(rst_om))

    n_jk <- nrow(data)
    dta_psbor_jk <- dta_psbor
    for (i_jk in 1:n_jk) {
        dta_psbor_jk$data <- data[-i_jk,]
        rst_jk <- get_ps_lrk_rmst(dta_psbor_jk, v_outcome = v_outcome,
                                  v_event = v_event, v_time = v_time,
                                  f_stratum = f_stratum,
                                  f_overall_est = f_overall_est, ...)
        sdf_om <- sdf_om + (rst_jk$Effect$Overall_Estimate$Mean - rst_om)^2
    }

    ## update rst
    nc_jk <- (n_jk - 1) / n_jk
    rst$Effect$Overall_Estimate$StdErr <- sqrt(sdf_om * nc_jk)

    ## return
    return(rst)
}


## Complex Jackknife

#' Get complex JK estimates for log-rank and RMST estimations between two arms
#' for RCT augmenting control
#'
#' @noRd
#'
get_ps_lrk_rmst_cjk <- function(dta_psbor,
                                v_outcome     = NULL,
                                v_event       = NULL,
                                v_time        = NULL,
                                f_stratum     = get_surv_stratum_lrk,
                                f_overall_est = get_overall_est_wostderr,
                                ...) {

    ## prepare data
    data_org <- dta_psbor$data
    Call_arg <- dta_psbor$Call_arg
    Call_fml <- dta_psbor$Call_fml
    for (i_call in 1:length(Call_arg)) {
        Call_arg[[i_call]][[".drop_arg_fml"]] <- TRUE
    }

    ## main estimates
    rst <- get_ps_lrk_rmst(dta_psbor, v_outcome = v_outcome,
                           v_event = v_event, v_time = v_time,
                           f_stratum = f_stratum,
                           f_overall_est = f_overall_est, ...)

    ## Complex JK stderr
    rst_om <- rst$Effect$Overall_Estimate$Mean
    sdf_om <- rep(0, length(rst_om))

    ## Complex JK
    n_jk <- nrow(data_org)
    n_jk_noerr <- n_jk

    for (i_jk in 1:n_jk) {
        rst_try <- try({
            ## get the jk dataset from the original data (without trimming).
            tmp_rst_jk <- data_org[-i_jk,]

            ## repeat all saved psrwe steps prior to get_ps_cl_km() call.
            ## > dta_ps <- psrwe_est(data, ...)
            ## > dta_ps_mat <- psrwe_match(dta_ps, ...)
            ## > dta_psbor_jk <- psrwe_match(dta_ps_mat, ...)
            for (i_call in 1:length(Call_fml)) {
                Call_arg[[i_call]][[1]] <- tmp_rst_jk
                tmp_rst_jk <- do.call(Call_fml[[i_call]], Call_arg[[i_call]])
            }
        }, silent = TRUE)

        if (inherits(rst_try, "try-error")) {
            n_jk_noerr <- n_jk_noerr - 1
            next
        }

        rst_jk <- get_ps_lrk_rmst(tmp_rst_jk, v_outcome = v_outcome,
                                  v_event = v_event, v_time = v_time,
                                  f_stratum = f_stratum,
                                  f_overall_est = f_overall_est, ...)
        sdf_om <- sdf_om + (rst_jk$Effect$Overall_Estimate$Mean - rst_om)^2
    }

    ## update rst
    nc_jk <- (n_jk_noerr - 1) / n_jk_noerr
    rst$Effect$Overall_Estimate$StdErr <- sqrt(sdf_om * nc_jk)

    ## return
    return(rst)
}

