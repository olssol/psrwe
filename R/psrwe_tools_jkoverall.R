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
                               f_stratum = f_stratum,
                               f_overall_est = f_overall_est, ...)

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

