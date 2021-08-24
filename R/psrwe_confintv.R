#' @title Confidence Interval for PS-Integrated Estimation
#'
#' Estimate the confidence interval for he PS-integrated approach.
#'
#' @param dta_psrst a returned object with class \code{RWE_PS_EST}
#' @param method_ci a method name for confidence interval (default Wald)
#' @param conf_type a type name of transformation for the confidence interal
#'        of PSKM approach (default log_log)
#' @param conf_int a two-sided level of confidence limits (default 0.95)
#' @param ... other options
#'
#' @return A list with class name \code{RWE_PS_EST}.
#'
#' @examples
#' \donttest{
#'
#' data(ex_dta)
#' dta_ps <- rwe_ps_est(ex_dta,
#'        v_covs = paste("V", 1:7, sep = ""),
#'        v_grp = "Group",
#'        cur_grp_level = "current")
#' ps_borrow <- rwe_ps_borrow(total_borrow = 30, dta_ps)
#' ps_rst <- rwe_ps_compl(ps_borrow, v_outcome = "Y_Con")
#' rst <- rwe_ps_ci(ps_rst)
#'}
#'
#' @export
#'
rwe_ps_ci <- function(dta_psrst,
                      method_ci = c("wald", "wilson"),
                      conf_type = c("log_log", "plain"),
                      conf_int = 0.95,
                      ...) {

    ## check
    stopifnot(inherits(dta_psrst,
                       what = get_rwe_class("PSRST")))

    stopifnot(dta_psrst$Method %in% c("ps_cl", "ps_km"))

    method_ci <- match.arg(method_ci)
    conf_type <- match.arg(conf_type)

    if (dta_psrst$Outcome_type != "binary") {
        method_ci <- "wald"
    }
    if (dta_psrst$Outcome_type != "tte") {
        conf_type <- NULL
    }

    ## prepare data
    is_rct <- dta_psrst$is_rct 

    n_ctl_s <- dta_psrst$Borrow$N_Current + dta_psrst$Borrow$N_Borrow
    n_ctl <- sum(n_ctl_s)

    if (is_rct) {
        n_eff_s <- n_ctl_s
        n_eff <- n_ctl
        n_ctl_s <- dta_psrst$Borrow$N_Cur_CTL + dta_psrst$Borrow$N_Borrow
        n_trt_s <- dta_psrst$Borrow$N_Cur_TRT
        n_ctl <- sum(n_ctl_s)
        n_trt <- sum(n_trt_s)
    }

    ## prepare for the return object
    rst_psci <- list(Control = NULL,
                     Treatment = NULL,
                     Effect = NULL,
                     Method = method_ci,
                     Conf_type = conf_type,
                     level_2sided = conf_int)

    ## by study type
    rst_psci$Control$Stratum_Estimate <-
        get_ci(dta_psrst$Control$Stratum_Estimate,
               method_ci, conf_type, conf_int, n_ctl_s, ...)
    rst_psci$Control$Overall_Estimate <-
        get_ci(dta_psrst$Control$Overall_Estimate,
               method_ci, conf_type, conf_int, n_ctl, ...)

    if (is_rct) {
        rst_psci$Treatment$Stratum_Estimate <-
            get_ci(dta_psrst$Treatment$Stratum_Estimate,
                   method_ci, conf_type, conf_int, n_trt_s, ...)
        rst_psci$Treatment$Overall_Estimate <-
            get_ci(dta_psrst$Treatment$Overall_Estimate,
                   method_ci, conf_type, conf_int, n_trt, ...)

        rst_psci$Effect$Stratum_Estimate <-
            get_ci(dta_psrst$Effect$Stratum_Estimate,
                   method_ci, conf_type, conf_int, n_eff_s, ...)
        rst_psci$Effect$Overall_Estimate <-
            get_ci(dta_psrst$Effect$Overall_Estimate,
                   method_ci, conf_type, conf_int, n_eff, ...)
    }

    ## return
    rst <- dta_psrst
    rst$CI <- rst_psci
    rst
}


#' @title Confidence interval (one arm)
#'
#' @noRd
get_ci <- function(x,
                   method_ci,
                   conf_type,
                   conf_int,
                   n,
                   ...) {

    if (method_ci == "wald") {
        if (is.null(conf_type) || is.na(conf_type)) {
            rst <- get_ci_wald(x$Mean, x$StdErr, conf_int, ...)
        } else {
            rst <- get_ci_km(x$Mean, x$StdErr, conf_int, conf_type, ...)
        }
    } else if (method_ci == "wilson") {
        rst <- get_ci_wilson(x$Mean, x$StdErr, conf_int, n, ...)
    } else {
        stop("Confidence interval method is not implemented.")
    }

    rst
}


#' @title Wilson Confidence interval for binary outcomes (one arm)
#'
#' @noRd
get_ci_wilson <- function(mean,
                          stderr,
                          conf_int = 0.95,
                          n,
                          method_stderr = c("original", "plain"),
                          ...) {

    z_alphad2 <- qnorm((1 - conf_int) / 2, lower.tail = FALSE)

    method_stderr <- match.arg(method_stderr)
    stderr <- switch(method_stderr,
                     original = {
                         sqrt(mean * (1 - mean) / n)
                     },
                     plain = {
                         stderr
                     })

    w_n <- n / (n + z_alphad2^2)
    w_z <- z_alphad2^2 / (n + z_alphad2^2)
    var_ws <- stderr^2 * w_n^2 + w_z^2 / (4 * z_alphad2^2)
    ci_ws_lb <- mean * w_n + w_z / 2 - z_alphad2 * sqrt(var_ws)
    ci_ws_ub <- mean * w_n + w_z / 2 + z_alphad2 * sqrt(var_ws)

    rst <- data.frame(Lower = ci_ws_lb, Upper = ci_ws_ub)
    rst
}


#' @title Wald Confidence interval for continous outcomes (one arm)
#'
#' @noRd
get_ci_wald <- function(mean,
                        stderr,
                        conf_int = 0.95,
                        ...) {

    z_alphad2 <- qnorm((1 - conf_int) / 2, lower.tail = FALSE)

    ci_wl_lb <- mean - stderr * z_alphad2 
    ci_wl_ub <- mean + stderr * z_alphad2 

    rst <- data.frame(Lower = ci_wl_lb, Upper = ci_wl_ub)
    rst
}


#' @title Wald Confidence interval for KM and time-to-event (one arm)
#'
#' @noRd
get_ci_km <- function(mean,
                      stderr,
                      conf_int = 0.95,
                      conf_type = c("log_log", "plain"),
                      ...) {

    conf_type <- match.arg(conf_type)
    z_alphad2 <- qnorm((1 - conf_int) / 2, lower.tail = FALSE)

    ci <- switch(conf_type,
                 log_log = {
                     log_S        <- log(mean)
                     se_log_log_S <- stderr / mean / log_S
                     A <- cbind(-z_alphad2 * se_log_log_S,
                                z_alphad2 * se_log_log_S)
                     ci <- mean^exp(A)
                 },
                 plain = cbind(mean - z_alphad2 * stderr,
                               mean + z_alphad2 * stderr)
                 )

    rst <- data.frame(Lower = ci[, 1], Upper = ci[, 2])
    rst
}

