#' @title Confidence/Credible Interval for PS-Integrated Estimation
#'
#' Estimate the confidence/credible interval for the PS-integrated approach.
#'
#' @param dta_psrst a returned object with class \code{PSRWE_EST}
#' @param method_ci a method name for confidence interval (default Wald)
#' @param conf_type a type name of transformation for the confidence interal
#'        of PSKM approach (default log_log)
#' @param conf_int a two-sided level of confidence/credible limits
#'        (default 0.95)
#' @param ... other options
#'
#' @return A list with class name \code{PSRWE_EST}.
#'
#' @details \code{method_ci = "wilson"} is for binary outcomes only.
#'
#' @examples
#' data(ex_dta)
#' dta_ps <- psrwe_est(ex_dta,
#'        v_covs = paste("V", 1:7, sep = ""),
#'        v_grp = "Group",
#'        cur_grp_level = "current")
#' ps_borrow <- psrwe_borrow(total_borrow = 30, dta_ps)
#' ps_rst <- psrwe_compl(ps_borrow, v_outcome = "Y_Con")
#' rst <- psrwe_ci(ps_rst)
#' rst
#'
#' @export
#'
psrwe_ci <- function(dta_psrst,
                     method_ci = c("wald", "wilson"),
                     conf_type = c("log_log", "plain"),
                     conf_int = 0.95,
                     ...) {

    ## check
    stopifnot(inherits(dta_psrst,
                       what = get_rwe_class("ANARST")))

    stopifnot(dta_psrst$Method %in% c("ps_pp", "ps_cl", "ps_km"))

    method_ci <- match.arg(method_ci)
    outcome_type <- dta_psrst$Outcome_type

    stopifnot(!(method_ci == "wilson" && outcome_type != "binary"))
    conf_type <- ifelse(outcome_type == "tte", match.arg(conf_type), NA)

    ## get ci by method
    if (dta_psrst$Method == "ps_pp") {
        rst_psci <- get_psci_bayesian(dta_psrst,
                                      conf_int,
                                      ...)
    } else {
        rst_psci <- get_psci_freq(dta_psrst,
                                  method_ci,
                                  conf_type,
                                  conf_int,
                                  ...)
    }


    ## return
    rst <- dta_psrst
    rst$CI <- rst_psci

    return(rst)
}




#' @title Bayesian credible interval
#'
#' @noRd
get_psci_bayesian <- function(dta_psrst,
                              conf_int,
                              ...) {

    ## prepare data
    is_rct <- dta_psrst$is_rct 

    ## prepare for the return object
    rst_psci <- list(Control = NULL,
                     Treatment = NULL,
                     Effect = NULL,
                     Method_ci = "credible interval",
                     Conf_type = NA,
                     Conf_int = conf_int)

    ## by study type
    rst_psci$Control$Stratum_Estimate <-
        get_bci(dta_psrst$Control$Stratum_Samples,
                conf_int = conf_int,
                ...)
    rst_psci$Control$Overall_Estimate <-
        get_bci(dta_psrst$Control$Overall_Samples,
                conf_int = conf_int,
                ...)

    if (is_rct) {
        rst_psci$Treatment$Stratum_Estimate <-
            get_bci(dta_psrst$Treatment$Stratum_Samples,
                    conf_int = conf_int,
                    ...)
        rst_psci$Treatment$Overall_Estimate <-
            get_bci(dta_psrst$Treatment$Overall_Samples,
                    conf_int = conf_int,
                    ...)

        rst_psci$Effect$Stratum_Estimate <-
            get_bci(dta_psrst$Effect$Stratum_Samples,
                    conf_int = conf_int,
                    ...)
        rst_psci$Effect$Overall_Estimate <-
            get_bci(dta_psrst$Effect$Overall_Samples,
                    conf_int = conf_int,
                    ...)
    }

    return(rst_psci)
}


#' @title Credible interval
#'
#' @noRd
get_bci <- function(x,
                    conf_int = 0.95,
                    ...) {

    q_alphad2 <- (1 - conf_int) / 2

    if (is.vector(x)) {
        bci_lb <- quantile(x, q_alphad2, names = FALSE)
        bci_ub <- quantile(x, 1 - q_alphad2, names = FALSE)
    } else {
        bci_lb <- apply(x, 1, quantile, q_alphad2) 
        bci_ub <- apply(x, 1, quantile, 1 - q_alphad2) 
    }

    rst <- data.frame(Lower = bci_lb, Upper = bci_ub)

    return(rst)
}




#' @title Frequentist confidence interval
#'
#' @noRd
get_psci_freq <- function(dta_psrst,
                          method_ci,
                          conf_type,
                          conf_int,
                          ...) {

    ## check
    outcome_type <- dta_psrst$Outcome_type

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
                     Method_ci = method_ci,
                     Conf_type = conf_type,
                     Conf_int = conf_int)

    ## by study type
    rst_psci$Control$Stratum_Estimate <-
        get_ci(dta_psrst$Control$Stratum_Estimate,
               n_ctl_s,
               method_ci = method_ci,
               conf_type = conf_type,
               conf_int = conf_int,
               ...)
    rst_psci$Control$Overall_Estimate <-
        get_ci(dta_psrst$Control$Overall_Estimate,
               n_ctl,
               method_ci = method_ci,
               conf_type = conf_type,
               conf_int = conf_int,
               ...)

    if (is_rct) {
        rst_psci$Treatment$Stratum_Estimate <-
            get_ci(dta_psrst$Treatment$Stratum_Estimate,
                   n_trt_s,
                   method_ci = method_ci,
                   conf_type = conf_type,
                   conf_int = conf_int,
                   ...)
        rst_psci$Treatment$Overall_Estimate <-
            get_ci(dta_psrst$Treatment$Overall_Estimate,
                   n_trt,
                   method_ci = method_ci,
                   conf_type = conf_type,
                   conf_int = conf_int,
                   ...)

        if (method_ci == "wald") {
            rst_psci$Effect$Stratum_Estimate <-
                get_ci(dta_psrst$Effect$Stratum_Estimate,
                       n_eff_s,
                       method_ci = method_ci,
                       conf_type = conf_type,
                       conf_int = conf_int,
                       ...)

            if (outcome_type != "tte" || conf_type == "plain") {
                rst_psci$Effect$Overall_Estimate <-
                    get_ci(dta_psrst$Effect$Overall_Estimate,
                           n_eff,
                           method_ci = method_ci,
                           conf_type = conf_type,
                           conf_int = conf_int,
                           ...)
            } else {
                ## TODO: Difference of two KMs with log-log.
                rst_psci$Effect$Overall_Estimate <- NA
            }
        } else {
            rst_psci$Effect$Stratum_Estimate <-
                get_ci_2arms(dta_psrst$Treatment$Stratum_Estimate,
                             dta_psrst$Control$Stratum_Estimate,
                             n_trt_s,
                             n_ctl_s,
                             method_ci = method_ci,
                             conf_type = conf_type,
                             conf_int = conf_int,
                             ...)

            ## TODO: Derive score method
            rst_psci$Effect$Overall_Estimate <- NA
        }
    }

    return(rst_psci)
}


#' @title Confidence interval (one arm)
#'
#' @noRd
get_ci <- function(x,
                   n = NULL,
                   method_ci = c("wald", "wilson"),
                   conf_type = c("log_log", "plain"),
                   conf_int = 0.95,
                   ...) {

    if (method_ci == "wald") {
        if (!any("T" %in% names(x))) {
            rst <- get_ci_wald(x$Mean,
                               x$StdErr,
                               conf_int = conf_int,
                               ...)
        } else {
            rst <- get_ci_km(x$Mean,
                             x$StdErr,
                             conf_type = conf_type,
                             conf_int = conf_int,
                             ...)
        }
    } else if (method_ci == "wilson") {
        rst <- get_ci_wilson(x$Mean,
                             x$StdErr,
                             n,
                             conf_int = conf_int,
                             ...)
    } else {
        stop("Confidence interval method is not implemented.")
    }

    return(rst)
}


#' @title Wald CI for continous outcomes (one arm)
#'
#' @noRd
get_ci_wald <- function(mean,
                        stderr,
                        conf_int = 0.95,
                        ...) {

    z_alphad2 <- qnorm((1 - conf_int) / 2, lower.tail = FALSE)

    ci_wl_lb <- mean - z_alphad2 * stderr
    ci_wl_ub <- mean + z_alphad2 * stderr

    rst <- data.frame(Lower = ci_wl_lb, Upper = ci_wl_ub)

    return(rst)
}


#' @title Wilson CI for binary outcomes (one arm)
#'
#' @noRd
get_ci_wilson <- function(mean,
                          stderr,
                          n,
                          conf_int = 0.95,
                          method_stderr = c("original", "plain"),
                          ...) {

    z_alphad2 <- qnorm((1 - conf_int) / 2, lower.tail = FALSE)

    method_stderr <- match.arg(method_stderr)
    if (method_stderr == "original") {
        stderr <- sqrt(mean * (1 - mean) / n)
    }

    two_a <- 2 * (n + z_alphad2^2)
    minus_b <- 2 * n * mean + z_alphad2^2
    b2_m_4ac <- (z_alphad2^2 + 4 * n * (stderr^2 * n)) * z_alphad2^2
    sqrt_b2_m_4ac <- sqrt(b2_m_4ac)
    ci_ws_lb <- (minus_b - sqrt_b2_m_4ac) / two_a
    ci_ws_ub <- (minus_b + sqrt_b2_m_4ac) / two_a

    rst <- data.frame(Lower = ci_ws_lb, Upper = ci_ws_ub)

    return(rst)
}


#' @title Wald CI for KM and time-to-event (one arm)
#'
#' @noRd
get_ci_km <- function(mean,
                      stderr,
                      conf_type = c("log_log", "plain"),
                      conf_int = 0.95,
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

    return(rst)
}


#' @title Confidence interval (two arms, independent)
#'
#' @noRd
get_ci_2arms <- function(x,
                         x_2,
                         n = NULL,
                         n_2 = NULL,
                         method_ci = c("wald", "wilson"),
                         conf_type = c("log_log", "plain"),
                         conf_int = 0.95,
                         ...) {

    if (method_ci == "wald") {
        if (!any("T" %in% names(x))) {
            rst <- get_ci_wald_2arms(x$Mean,
                                     x$StdErr,
                                     x_2$Mean,
                                     x_2$StdErr,
                                     conf_int = conf_int,
                                     ...)
        } else {
            rst <- get_ci_km_2arms(x$Mean,
                                   x$StdErr,
                                   x_2$Mean,
                                   x_2$StdErr,
                                   conf_int = conf_int,
                                   conf_type = conf_type,
                                   ...)
        }
    } else if (method_ci == "wilson") {
        rst <- get_ci_wilson_2arms(x$Mean,
                                   x$StdErr,
                                   x_2$Mean,
                                   x_2$StdErr,
                                   n,
                                   n_2,
                                   conf_int = conf_int,
                                   ...)
    } else {
        stop("Confidence interval method is not implemented.")
    }

    return(rst)
}


#' @title Wald CI for continous outcomes (two arms, independent)
#'
#' @noRd
get_ci_wald_2arms <- function(mean,
                              stderr,
                              mean_2,
                              stderr_2,
                              conf_int = 0.95,
                              ...) {

    z_alphad2 <- qnorm((1 - conf_int) / 2, lower.tail = FALSE)

    stderr_p <- sqrt(stderr^2 + stderr_2^2)
    ci_wl_lb <- (mean - mean_2) - z_alphad2 * stderr_p
    ci_wl_ub <- (mean - mean_2) + z_alphad2 * stderr_p

    rst <- data.frame(Lower = ci_wl_lb, Upper = ci_wl_ub)

    return(rst)
}


#' @title Wilson CI for binary outcomes (two arms, independent)
#'
#' @noRd
get_ci_wilson_2arms <- function(mean,
                                stderr,
                                mean_2,
                                stderr_2,
                                n,
                                n_2,
                                conf_int = 0.95,
                                method_stderr = c("original", "plain"),
                                ...) {

    z_alphad2 <- qnorm((1 - conf_int) / 2, lower.tail = FALSE)

    method_stderr <- match.arg(method_stderr)
    if (method_stderr == "original") {
        stderr <- sqrt(mean * (1 - mean) / n)
        stderr_2 <- sqrt(mean_2 * (1 - mean_2) / n_2)
    }

    get_root <- function(p, se, n, z) {
        two_a <- 2 * (n + z^2)
        minus_b <- 2 * n * p + z^2
        b2_m_4ac <- (z^2 + 4 * n * (se^2 * n)) * z^2
        sqrt_b2_m_4ac <- sqrt(b2_m_4ac)
        lb <- (minus_b - sqrt_b2_m_4ac) / two_a
        ub <- (minus_b + sqrt_b2_m_4ac) / two_a
        cbind(lb, ub)
    }

    rt <- get_root(mean, stderr, n, z_alphad2)
    rt_2 <- get_root(mean_2, stderr_2, n_2, z_alphad2)

    stderr_d <- sqrt(rt[, 1] * (1 - rt[, 1]) / n +
                     rt_2[, 2] * (1 - rt_2[, 2]) / n_2)
    stderr_e <- sqrt(rt[, 2] * (1 - rt[, 2]) / n +
                     rt_2[, 1] * (1 - rt_2[, 1]) / n_2)

    ci_ws_lb <- (mean - mean_2) - z_alphad2 * stderr_d
    ci_ws_ub <- (mean - mean_2) + z_alphad2 * stderr_e

    rst <- data.frame(Lower = ci_ws_lb, Upper = ci_ws_ub)

    return(rst)
}


#' @title Wald CI for KM and time-to-event (two arms, independent)
#'
#' @noRd
get_ci_km_2arms <- function(mean,
                            stderr,
                            mean_2,
                            stderr_2,
                            conf_int = 0.95,
                            conf_type = c("log_log", "plain"),
                            ...) {

    conf_type <- match.arg(conf_type)
    z_alphad2 <- qnorm((1 - conf_int) / 2, lower.tail = FALSE)

    ci <- switch(conf_type,
                 log_log = {
                     ## TODO: Difference of two KMs with log-log.
                     cbind(rep(NA, length(mean)),
                           rep(NA, length(mean)))
                 },
                 plain = {
                     stderr_p <- sqrt(stderr^2 + stderr_2^2)
                     cbind((mean - mean_2) - z_alphad2 * stderr_p,
                           (mean - mean_2) + z_alphad2 * stderr_p)
                 })

    rst <- data.frame(Lower = ci[, 1], Upper = ci[, 2])

    return(rst)
}

