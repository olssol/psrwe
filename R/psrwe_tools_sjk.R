#' Get simple JKoverall estimates for composite likelihood and survival
#'
#' @noRd
#'
get_ps_cl_km_sjk <- function(dta_psbor,
                             v_outcome     = NULL,
                             v_event       = NULL,
                             v_time        = NULL,
                             f_stratum     = get_cl_stratum,
                             f_overall_est = get_overall_est_wostderr,
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

    ## simple JK stderr
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

    n_jk <- nrow(data)
    dta_psbor_jk <- dta_psbor
    for (i_jk in 1:n_jk) {
        dta_psbor_jk$data <- data[-i_jk,]
        rst_jk <- get_ps_cl_km(dta_psbor_jk, v_outcome = v_outcome,
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
    nc_jk <- (n_jk - 1) / n_jk
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


#' Get simple Bootstrap estimates for composite likelihood and survival
#'
#' @noRd
#'
get_ps_cl_km_sbs <- function(dta_psbor,
                             v_outcome     = NULL,
                             v_event       = NULL,
                             v_time        = NULL,
                             f_stratum     = get_cl_stratum,
                             f_overall_est = get_overall_est_wostderr,
                             n_bootstrap   = 200,
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

    ## simple Bootstrap stderr
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

        rst_bs <- get_ps_cl_km(dta_psbor_bs, v_outcome = v_outcome,
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
    nc_bootstrap <- 1 / n_bootstrap
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


