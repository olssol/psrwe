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
    Call_fml <- dta_psbor$Call_fml

    ## main estimates
    rst <- get_ps_cl_km(dta_psbor, v_outcome = v_outcome,
                        v_event = v_event, v_time = v_time,
                        f_stratum = f_stratum,
                        f_overall_est = f_overall_est, ...)

    ## Complex JK stderr
    rstom_ctl <- rst$Control$Overall_Estimate$Mean
    sdf_ctl <- rep(0, length(rstom_ctl))

    if (is_rct) {
        rstom_trt <- rst$Treatment$Overall_Estimate$Mean
        rstom_eff <- rst$Effect$Overall_Estimate$Mean
        sdf_trt <- rep(0, length(rstom_trt))
        sdf_eff <- rep(0, length(rstom_eff))
    }

    ## Complex JK
    n_jk <- nrow(data_org)
    n_jk_noerr <- n_jk

    for (i_jk in 1:n_jk) {
        ## Avoid interruptions by errors during PS steps.
        rst_try <- try({
            ## Get the jk dataset from the original data (without trimming).
            data_jk <- data_org[-i_jk,]

            ## Repeat all saved psrwe steps prior to get_ps_cl_km() call.
            ## dta_ps <- psrwe_est(data, ...)
            ## dta_ps_mat <- psrwe_match(dta_ps, ...)
            ## dta_psbor_jk <- psrwe_match(dta_ps_mat, ...)

            name_rst <- as.character(Call_fml[[1]][[2]])
            assign(name_rst, data_jk)

            for (i_call in 1:(length(Call_fml) - 1)) {
              tmp_rst  <- eval(Call_fml[[i_call]])
              name_rst <- as.character(Call_fml[[i_call + 1]][[2]])
              assign(name_rst, tmp_rst)
            }

            dta_psbor_jk <- eval(Call_fml[[i_call + 1]])
        }, silent = TRUE)

        if (inherits(rst_try, "try-error")) {
            n_jk_noerr <- n_jk_noerr - 1
            next
        }

        ## Do the same as "jkoverall" method on the new "dta_psbor_jk".
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
    nc_jk <- (n_jk_noerr - 1) / n_jk_noerr
    rst$Control$Overall_Estimate$StdErr <- sqrt(sdf_ctl * nc_jk)

    if (is_rct) {
        rst$Treatment$Overall_Estimate$StdErr <- sqrt(sdf_trt * nc_jk)
        rst$Effect$Overall_Estimate$StdErr <- sqrt(sdf_eff * nc_jk)
    }

    ## return
    return(rst)
}

