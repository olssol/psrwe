#' Get complex JK estimates for composite likelihood and survival
#'
#' @noRd
#'
get_ps_cl_km_cjk <- function(dta_psbor,
                             v_outcome     = NULL,
                             v_event       = NULL,
                             v_time        = NULL,
                             f_stratum     = get_cl_stratum,
                             f_overall_est = get_overall_est_wostderr,
                             ...) {

    ## prepare data
    is_rct   <- dta_psbor$is_rct
    data_org <- dta_psbor$data
    Call_arg <- dta_psbor$Call_arg
    Call_fml <- dta_psbor$Call_fml
    for (i_call in 1:length(Call_arg)) {
        Call_arg[[i_call]][[".drop_arg_fml"]] <- TRUE
    }

    ## main estimates
    rst <- get_ps_cl_km(dta_psbor, v_outcome = v_outcome,
                        v_event = v_event, v_time = v_time,
                        f_stratum = f_stratum,
                        f_overall_est = f_overall_est, ...)

    ## complex JK stderr
    rst_sm_ctl <- rst$Control$Stratum_Estimate$Mean
    rst_om_ctl <- rst$Control$Overall_Estimate$Mean
    sdf_sm_ctl <- rep(0, length(rst_sm_ctl))
    sdf_om_ctl <- rep(0, length(rst_om_ctl))

    if (is_rct) {
        rst_sm_trt <- rst$Treatment$Stratum_Estimate$Mean
        rst_om_trt <- rst$Treatment$Overall_Estimate$Mean
        sdf_sm_trt <- rep(0, length(rst_sm_trt))
        sdf_om_trt <- rep(0, length(rst_om_trt))

        rst_sm_eff <- rst$Effect$Stratum_Estimate$Mean
        rst_om_eff <- rst$Effect$Overall_Estimate$Mean
        sdf_sm_eff <- rep(0, length(rst_sm_eff))
        sdf_om_eff <- rep(0, length(rst_om_eff))
    }

    ## complex JK
    n_jk <- nrow(data_org)
    n_jk_noerr <- n_jk

    for (i_jk in 1:n_jk) {
        rst_try <- try({
            ## get the jk dataset from the original data (without trimming).
            tmp_rst_jk <- data_org[-i_jk,]

            ## repeat all saved psrwe steps prior to get_ps_cl_km() call.
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
        rst_jk <- get_ps_cl_km(tmp_rst_jk, v_outcome = v_outcome,
                               v_event = v_event, v_time = v_time,
                               f_stratum = f_stratum,
                               f_overall_est = f_overall_est, ...)

        sdf_sm_ctl <- sdf_sm_ctl + (rst_jk$Control$Stratum_Estimate$Mean -
                                    rst_sm_ctl)^2
        sdf_om_ctl <- sdf_om_ctl + (rst_jk$Control$Overall_Estimate$Mean -
                                    rst_om_ctl)^2

        if (is_rct) {
            sdf_sm_trt <- sdf_sm_trt + (rst_jk$Treatment$Stratum_Estimate$Mean -
                                        rst_sm_trt)^2
            sdf_om_trt <- sdf_om_trt + (rst_jk$Treatment$Overall_Estimate$Mean -
                                        rst_om_trt)^2

            sdf_sm_eff <- sdf_sm_eff + (rst_jk$Effect$Stratum_Estimate$Mean -
                                        rst_sm_eff)^2
            sdf_om_eff <- sdf_om_eff + (rst_jk$Effect$Overall_Estimate$Mean -
                                        rst_om_eff)^2
        }
    }

    ## update rst
    stopifnot(n_jk_noerr > 1)
    nc_jk <- (n_jk_noerr - 1) / n_jk_noerr
    rst$Control$Stratum_Estimate$StdErr <- sqrt(sdf_sm_ctl * nc_jk)
    rst$Control$Overall_Estimate$StdErr <- sqrt(sdf_om_ctl * nc_jk)

    if (is_rct) {
        rst$Treatment$Stratum_Estimate$StdErr <- sqrt(sdf_sm_trt * nc_jk)
        rst$Treatment$Overall_Estimate$StdErr <- sqrt(sdf_om_trt * nc_jk)

        rst$Effect$Stratum_Estimate$StdErr <- sqrt(sdf_sm_eff * nc_jk)
        rst$Effect$Overall_Estimate$StdErr <- sqrt(sdf_om_eff * nc_jk)
    }

    ## return
    return(rst)
}


#' Get complex Bootstrap estimates for composite likelihood and survival
#'
#' @noRd
#'
get_ps_cl_km_cbs <- function(dta_psbor,
                             v_outcome     = NULL,
                             v_event       = NULL,
                             v_time        = NULL,
                             f_stratum     = get_cl_stratum,
                             f_overall_est = get_overall_est_wostderr,
                             n_bootstrap   = 200,
                             ...) {

    ## prepare data
    is_rct   <- dta_psbor$is_rct
    data_org <- dta_psbor$data
    data_org <- data_org[!is.na(data_org[["_strata_"]]), ]
    Call_arg <- dta_psbor$Call_arg
    Call_fml <- dta_psbor$Call_fml
    for (i_call in 1:length(Call_arg)) {
        Call_arg[[i_call]][[".drop_arg_fml"]] <- TRUE
    }

    ## main estimates
    rst <- get_ps_cl_km(dta_psbor, v_outcome = v_outcome,
                        v_event = v_event, v_time = v_time,
                        f_stratum = f_stratum,
                        f_overall_est = f_overall_est, ...)

    ## complex Bootstrap stderr
    rst_sm_ctl <- rst$Control$Stratum_Estimate$Mean
    rst_om_ctl <- rst$Control$Overall_Estimate$Mean
    sdf_sm_ctl <- rep(0, length(rst_sm_ctl))
    sdf_om_ctl <- rep(0, length(rst_om_ctl))

    if (is_rct) {
        rst_sm_trt <- rst$Treatment$Stratum_Estimate$Mean
        rst_om_trt <- rst$Treatment$Overall_Estimate$Mean
        sdf_sm_trt <- rep(0, length(rst_sm_trt))
        sdf_om_trt <- rep(0, length(rst_om_trt))

        rst_sm_eff <- rst$Effect$Stratum_Estimate$Mean
        rst_om_eff <- rst$Effect$Overall_Estimate$Mean
        sdf_sm_eff <- rep(0, length(rst_sm_eff))
        sdf_om_eff <- rep(0, length(rst_om_eff))
    }

    ## get id by unique stratum
    ustrata <- unique(data_org[, c("_grp_", "_arm_", "_strata_")])
    ustrata <- ustrata[!is.na(ustrata$"_strata"),]
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

            ## repeat all saved psrwe steps prior to get_ps_cl_km() call.
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
        rst_bs <- get_ps_cl_km(tmp_rst_bs, v_outcome = v_outcome,
                               v_event = v_event, v_time = v_time,
                               f_stratum = f_stratum,
                               f_overall_est = f_overall_est, ...)

        sdf_sm_ctl <- sdf_sm_ctl + (rst_bs$Control$Stratum_Estimate$Mean -
                                    rst_sm_ctl)^2
        sdf_om_ctl <- sdf_om_ctl + (rst_bs$Control$Overall_Estimate$Mean -
                                    rst_om_ctl)^2

        if (is_rct) {
            sdf_sm_trt <- sdf_sm_trt + (rst_bs$Treatment$Stratum_Estimate$Mean -
                                        rst_sm_trt)^2
            sdf_om_trt <- sdf_om_trt + (rst_bs$Treatment$Overall_Estimate$Mean -
                                        rst_om_trt)^2

            sdf_sm_eff <- sdf_sm_eff + (rst_bs$Effect$Stratum_Estimate$Mean -
                                        rst_sm_eff)^2
            sdf_om_eff <- sdf_om_eff + (rst_bs$Effect$Overall_Estimate$Mean -
                                        rst_om_eff)^2
        }
    }

    ## update rst
    stopifnot(n_bootstrap_noerr > 1)
    nc_bootstrap <- 1 / n_bootstrap_noerr
    rst$Control$Stratum_Estimate$StdErr <- sqrt(sdf_sm_ctl * nc_bootstrap)
    rst$Control$Overall_Estimate$StdErr <- sqrt(sdf_om_ctl * nc_bootstrap)

    if (is_rct) {
        rst$Treatment$Stratum_Estimate$StdErr <- sqrt(sdf_sm_trt * nc_bootstrap)
        rst$Treatment$Overall_Estimate$StdErr <- sqrt(sdf_om_trt * nc_bootstrap)

        rst$Effect$Stratum_Estimate$StdErr <- sqrt(sdf_sm_eff * nc_bootstrap)
        rst$Effect$Overall_Estimate$StdErr <- sqrt(sdf_om_eff * nc_bootstrap)
    }

    ## return
    return(rst)
}


