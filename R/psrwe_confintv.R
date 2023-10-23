#' @title Confidence/Credible Interval for PS-Integrated Estimation
#'
#' @description
#' Estimate the confidence/credible interval for the PS-integrated approach.
#'
#' @param dta_psrst A returned object with class \code{PSRWE_EST}
#' @param method_ci A method name for confidence interval (default wald)
#' @param conf_int A two-sided level of confidence/credible limits
#'        (default 0.95)
#' @param conf_type A type name of transformation for the confidence interval
#'        of PSKM approach
#' @param ... Other options
#'
#' @return A list with class name \code{PSRWE_EST}.
#'
#' @details \code{method_ci = "wilson"} is for binary outcomes only.
#'          \code{conf_type = "log_log"} is for \code{ps_km} only.
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
                     conf_int = 0.95,
                     conf_type = c("log_log", "plain"),
                     ...) {

    ## check
    stopifnot(inherits(dta_psrst,
                       what = get_rwe_class("ANARST")))

    stopifnot(dta_psrst$Method %in% get_rwe_class("ANAMETHOD"))

    method_ci <- match.arg(method_ci)
    outcome_type <- dta_psrst$Outcome_type

    stopifnot(!(method_ci == "wilson" && outcome_type != "binary"))
    conf_type <- ifelse(outcome_type == "tte", match.arg(conf_type), NA)

    ## get ci by method
    if (dta_psrst$Method == "ps_pp") {
        rst_psci <- get_psci_bayesian(dta_psrst,
                                      conf_int,
                                      ...)
    } else if (dta_psrst$Method == "ps_cl") {
        rst_psci <- get_psci_freq(dta_psrst,
                                  method_ci,
                                  conf_int,
                                  ...)
    } else if (dta_psrst$Method %in% get_rwe_class("ANAMETHOD_KM")) {
        if (dta_psrst$Method %in% c("ps_lrk", "ps_rmst")) {
             conf_type <- "plain"
        }
        rst_psci <- get_psci_km(dta_psrst,
                                conf_int,
                                conf_type,
                                ...)
    }

    ## return
    rst <- dta_psrst
    rst$CI <- rst_psci

    return(rst)
}




#' @title Bayesian credible interval (PSPP)
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
                     Conf_int = conf_int,
                     Conf_stderr = NA)

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




#' @title Frequentist confidence interval (PSCL)
#'
#' @noRd
get_psci_freq <- function(dta_psrst,
                          method_ci,
                          conf_int,
                          ...) {

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
                     Conf_type = NA,
                     Conf_int = conf_int,
                     Conf_stderr = NA)

    if (method_ci == "wilson") {
        tmp_conf_stderr <- list(...)$conf_stderr
        if (is.null(tmp_conf_stderr)){
            tmp_conf_stderr <- "original"
        }
        rst_psci$Conf_stderr <- tmp_conf_stderr
    }

    ## by study type
    rst_psci$Control$Stratum_Estimate <-
        get_fci(dta_psrst$Control$Stratum_Estimate,
                n_ctl_s,
                method_ci = method_ci,
                conf_int = conf_int,
                ...)
    rst_psci$Control$Overall_Estimate <-
        get_fci(dta_psrst$Control$Overall_Estimate,
                n_ctl,
                method_ci = method_ci,
                conf_int = conf_int,
                ...)

    if (is_rct) {
        rst_psci$Treatment$Stratum_Estimate <-
            get_fci(dta_psrst$Treatment$Stratum_Estimate,
                    n_trt_s,
                    method_ci = method_ci,
                    conf_int = conf_int,
                    ...)
        rst_psci$Treatment$Overall_Estimate <-
            get_fci(dta_psrst$Treatment$Overall_Estimate,
                    n_trt,
                    method_ci = method_ci,
                    conf_int = conf_int,
                    ...)

        if (method_ci == "wald") {
            rst_psci$Effect$Stratum_Estimate <-
                get_fci(dta_psrst$Effect$Stratum_Estimate,
                        n_eff_s,
                        method_ci = method_ci,
                        conf_int = conf_int,
                        ...)
            rst_psci$Effect$Overall_Estimate <-
                get_fci(dta_psrst$Effect$Overall_Estimate,
                        n_eff,
                        method_ci = method_ci,
                        conf_int = conf_int,
                        ...)
        } else {
            rst_psci$Effect$Stratum_Estimate <-
                get_fci_2arms(dta_psrst$Treatment$Stratum_Estimate,
                              dta_psrst$Control$Stratum_Estimate,
                              n_trt_s,
                              n_ctl_s,
                              method_ci = method_ci,
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
get_fci <- function(x,
                    n = NULL,
                    method_ci = c("wald", "wilson"),
                    conf_int = 0.95,
                    ...) {

    method_ci <- match.arg(method_ci)

    if (method_ci == "wald") {
        rst <- get_fci_wald(x$Mean,
                            x$StdErr,
                            conf_int = conf_int,
                            ...)
    } else if (method_ci == "wilson") {
        rst <- get_fci_wilson(x$Mean,
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
get_fci_wald <- function(mean,
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
get_fci_wilson <- function(mean,
                           stderr,
                           n,
                           conf_int = 0.95,
                           conf_stderr = c("original", "plain"),
                           ...) {

    z_alphad2 <- qnorm((1 - conf_int) / 2, lower.tail = FALSE)

    conf_stderr <- match.arg(conf_stderr)
    if (conf_stderr == "original") {
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




#' @title Confidence interval (two arms, independent)
#'
#' @noRd
get_fci_2arms <- function(x,
                          x_2,
                          n = NULL,
                          n_2 = NULL,
                          method_ci = c("wald", "wilson"),
                          conf_int = 0.95,
                          ...) {

    method_ci <- match.arg(method_ci)

    if (method_ci == "wald") {
        rst <- get_fci_2arms_wald(x$Mean,
                                  x$StdErr,
                                  x_2$Mean,
                                  x_2$StdErr,
                                  conf_int = conf_int,
                                  ...)
    } else if (method_ci == "wilson") {
        rst <- get_fci_2arms_wilson(x$Mean,
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
get_fci_2arms_wald <- function(mean,
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
get_fci_2arms_wilson <- function(mean,
                                 stderr,
                                 mean_2,
                                 stderr_2,
                                 n,
                                 n_2,
                                 conf_int = 0.95,
                                 conf_stderr = c("original", "plain"),
                                 ...) {

    z_alphad2 <- qnorm((1 - conf_int) / 2, lower.tail = FALSE)

    conf_stderr <- match.arg(conf_stderr)
    if (conf_stderr == "original") {
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




#' @title KM confidence interval (PSKM)
#'
#' @noRd
get_psci_km <- function(dta_psrst,
                        conf_int,
                        conf_type,
                        ...) {

    ## prepare data
    is_rct <- dta_psrst$is_rct 

    ## prepare for the return object
    rst_psci <- list(Control = NULL,
                     Treatment = NULL,
                     Effect = NULL,
                     Method_ci = "wald",
                     Conf_type = conf_type,
                     Conf_int = conf_int,
                     Conf_stderr = NA)

    ## by study type
    if (exists("Control", dta_psrst) && !is.null(dta_psrst$Control)) {
        rst_psci$Control$Stratum_Estimate <-
            get_kmci(dta_psrst$Control$Stratum_Estimate,
                     conf_int = conf_int,
                     conf_type = conf_type,
                     ...)
        rst_psci$Control$Overall_Estimate <-
            get_kmci(dta_psrst$Control$Overall_Estimate,
                     conf_int = conf_int,
                     conf_type = conf_type,
                     ...)
        rst_psci$Control$Conf_type <- conf_type
    }

    if (exists("Treatment", dta_psrst) && !is.null(dta_psrst$Treatment)) {
        rst_psci$Treatment$Stratum_Estimate <-
            get_kmci(dta_psrst$Treatment$Stratum_Estimate,
                     conf_int = conf_int,
                     conf_type = conf_type,
                     ...)
        rst_psci$Treatment$Overall_Estimate <-
            get_kmci(dta_psrst$Treatment$Overall_Estimate,
                     conf_int = conf_int,
                     conf_type = conf_type,
                     ...)
        rst_psci$Treatment$Conf_type <- conf_type
    }

    if (exists("Effect", dta_psrst) && !is.null(dta_psrst$Effect)) {
        rst_psci$Effect$Stratum_Estimate <-
            get_kmci(dta_psrst$Effect$Stratum_Estimate,
                     conf_int = conf_int,
                     conf_type = "plain",
                     ...)
        rst_psci$Effect$Overall_Estimate <-
            get_kmci(dta_psrst$Effect$Overall_Estimate,
                     conf_int = conf_int,
                     conf_type = "plain",
                     ...)
        rst_psci$Effect$Conf_type <- "plain"

        ## Only "plain" ci is available for treatment effect in RCT.
        rst_psci$Conf_type <- "plain"
    }

    return(rst_psci)
}


#' @title KM confidence interval (one arm)
#'
#' @noRd
get_kmci <- function(x,
                     conf_int = 0.95,
                     conf_type = c("log_log", "plain"),
                     ...) {

    conf_type <- match.arg(conf_type)

    if (conf_type == "log_log" || conf_type == "plain") {
        rst <- get_kmci_wald(x$Mean,
                             x$StdErr,
                             conf_int = conf_int,
                             conf_type = conf_type,
                             ...)
    } else {
        stop("Confidence interval type is not implemented.")
    }

    return(rst)
}


#' @title Wald CI for KM and time-to-event (one arm)
#'
#' @noRd
get_kmci_wald <- function(mean,
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

    return(rst)
}

