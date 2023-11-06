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

#' Get simple JKoverall estimates for log-rank and RMST estimations between two arms
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

    ## complex JK stderr
    rst_om <- rst$Effect$Overall_Estimate$Mean
    sdf_om <- rep(0, length(rst_om))

    ## complex JK
    n_jk <- nrow(data_org)
    n_jk_noerr <- n_jk

    for (i_jk in 1:n_jk) {
        rst_try <- try({
            ## get the jk dataset from the original data (without trimming).
            tmp_rst_jk <- data_org[-i_jk,]

            ## repeat all saved psrwe steps prior to get_ps_lrk_rmst() call.
            ## > dta_ps <- psrwe_est(data, ...)
            ## > dta_ps_mat <- psrwe_match(dta_ps, ...)
            ## > dta_psbor_jk <- psrwe_borrow(dta_ps_mat, ...)
            for (i_call in 1:length(Call_fml)) {
                Call_arg[[i_call]][[1]] <- tmp_rst_jk
                tmp_rst_jk <- do.call(Call_fml[[i_call]], Call_arg[[i_call]])
            }
        }, silent = TRUE)

        if (inherits(rst_try, "try-error")) {
            n_jk_noerr <- n_jk_noerr - 1
            next
        }

        ## do the same as "sjk" method on the new "tmp_rst_jk".
        rst_jk <- get_ps_lrk_rmst(tmp_rst_jk, v_outcome = v_outcome,
                                  v_event = v_event, v_time = v_time,
                                  f_stratum = f_stratum,
                                  f_overall_est = f_overall_est, ...)
        sdf_om <- sdf_om + (rst_jk$Effect$Overall_Estimate$Mean - rst_om)^2
    }

    ## update rst
    stopifnot(n_jk_noerr > 1)
    nc_jk <- (n_jk_noerr - 1) / n_jk_noerr
    rst$Effect$Overall_Estimate$StdErr <- sqrt(sdf_om * nc_jk)

    ## return
    return(rst)
}


## Simple Bootstrap

#' Get simple Bootstrap estimates for log-rank and RMST estimations between two arms
#' for RCT augmenting control
#'
#' @noRd
#'
get_ps_lrk_rmst_sbs <- function(dta_psbor,
                                v_outcome     = NULL,
                                v_event       = NULL,
                                v_time        = NULL,
                                f_stratum     = get_surv_stratum_lrk,
                                f_overall_est = get_overall_est_wostderr,
                                n_bootstrap   = 200,
                                ...) {

    ## prepare data
    data    <- dta_psbor$data
    data    <- data[!is.na(data[["_strata_"]]), ]

    ## main estimates
    rst <- get_ps_lrk_rmst(dta_psbor, v_outcome = v_outcome,
                           v_event = v_event, v_time = v_time,
                           f_stratum = f_stratum,
                           f_overall_est = f_overall_est, ...)

    ## bootstrap overall stderr
    rst_om <- rst$Effect$Overall_Estimate$Mean
    sdf_om <- rep(0, length(rst_om))


    ## get id by unique stratum
    ustrata <- unique(data[, c("_grp_", "_arm_", "_strata_")])
    n_ustrata <- nrow(ustrata)
    rst_id <- list()
    for (i_ustrata in 1:n_ustrata) {
        rst_id[[i_ustrata]] <-
            which(data$"_grp_" == ustrata[i_ustrata, "_grp_"] &
                  data$"_arm_" == ustrata[i_ustrata, "_arm_"] &
                  data$"_strata_" == ustrata[i_ustrata, "_strata_"])
    }

    dta_psbor_bs <- dta_psbor
    for (i_bs in 1:n_bootstrap) {
        ## sample by each unique stratum with replacement
        tmp_data <- data
        for (i_ustrata in 1:n_ustrata) {
            new_id <- sample(rst_id[[i_ustrata]], replace = TRUE)
            tmp_data[rst_id[[i_ustrata]],] <- data[new_id,]
        }
        dta_psbor_bs$data <- tmp_data

        rst_bs <- get_ps_lrk_rmst(dta_psbor_bs, v_outcome = v_outcome,
                                  v_event = v_event, v_time = v_time,
                                  f_stratum = f_stratum,
                                  f_overall_est = f_overall_est, ...)
        sdf_om <- sdf_om + (rst_bs$Effect$Overall_Estimate$Mean - rst_om)^2
    }

    ## update rst
    nc_bootstrap <- 1 / n_bootstrap
    rst$Effect$Overall_Estimate$StdErr <- sqrt(sdf_om * nc_bootstrap)

    ## return
    return(rst)
}


## Complex Bootstrap

#' Get complex Bootstrap estimates for log-rank and RMST estimations between two arms
#' for RCT augmenting control
#'
#' @noRd
#'
get_ps_lrk_rmst_cbs <- function(dta_psbor,
                                v_outcome     = NULL,
                                v_event       = NULL,
                                v_time        = NULL,
                                f_stratum     = get_surv_stratum_lrk,
                                f_overall_est = get_overall_est_wostderr,
                                n_bootstrap   = 200,
                                ...) {

    ## prepare data
    data_org <- dta_psbor$data
    data_org <- data_org[!is.na(data_org[["_strata_"]]), ]
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

    ## complex Bootstrap stderr
    rst_om <- rst$Effect$Overall_Estimate$Mean
    sdf_om <- rep(0, length(rst_om))

    ## get id by unique stratum
    ustrata <- unique(data_org[, c("_grp_", "_arm_", "_strata_")])
    n_ustrata <- nrow(ustrata)
    rst_id <- list()
    for (i_ustrata in 1:n_ustrata) {
        rst_id[[i_ustrata]] <-
            which(data_org$"_grp_" == ustrata[i_ustrata, "_grp_"] &
                  data_org$"_arm_" == ustrata[i_ustrata, "_arm_"] &
                  data_org$"_strata_" == ustrata[i_ustrata, "_strata_"])
    }

    ## complex Bootstrap
    n_bootstrap_noerr <- n_bootstrap
    for (i_bs in 1:n_bootstrap) {
        rst_try <- try({
            ## sample by each unique stratum with replacement
            tmp_rst_bs <- data_org
            for (i_ustrata in 1:n_ustrata) {
                new_id <- sample(rst_id[[i_ustrata]], replace = TRUE)
                tmp_rst_bs[rst_id[[i_ustrata]],] <- data_org[new_id,]
            }

            ## repeat all saved psrwe steps prior to get_ps_lrk_rmst() call.
            ## > dta_ps <- psrwe_est(data, ...)
            ## > dta_ps_mat <- psrwe_match(dta_ps, ...)
            ## > dta_psbor_bs <- psrwe_borrow(dta_ps_mat, ...)
            for (i_call in 1:length(Call_fml)) {
                Call_arg[[i_call]][[1]] <- tmp_rst_bs
                tmp_rst_bs <- do.call(Call_fml[[i_call]], Call_arg[[i_call]])
            }
        }, silent = TRUE)

        if (inherits(rst_try, "try-error")) {
            n_bootstrap_noerr <- n_bootstrap_noerr - 1
            next
        }

        ## do the same as "sbs" method on the new "tmp_rst_bs".
        rst_bs <- get_ps_lrk_rmst(tmp_rst_bs, v_outcome = v_outcome,
                                  v_event = v_event, v_time = v_time,
                                  f_stratum = f_stratum,
                                  f_overall_est = f_overall_est, ...)
        sdf_om <- sdf_om + (rst_bs$Effect$Overall_Estimate$Mean - rst_om)^2
    }

    ## update rst
    stopifnot(n_bootstrap_noerr > 1)
    nc_bootstrap <- 1 / n_bootstrap_noerr
    rst$Effect$Overall_Estimate$StdErr <- sqrt(sdf_om * nc_bootstrap)

    ## return
    return(rst)
}


