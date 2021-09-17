#' PS-Integrated Composite Likelihood Estimation
#'
#' Estimate the mean of the outcome based on PS-integrated composite likelihood
#' approach. Variance is estimated by Jack-Knife method. Applies to the case
#' when there is only one external data source.
#'
#' @inheritParams psrwe_powerp
#'
#' @param ... Parameters for \code{rwe_cl}
#'
#' @return A data frame with class name \code{PSRWE_RST}. It contains the
#'     composite estimation of the mean for each stratum as well as the
#'     jackknife estimation for each subject. The results should be further
#'     summarized by its S3 method \code{summary}.
#'
#' @examples
#' data(ex_dta)
#' dta_ps <- psrwe_est(ex_dta,
#'        v_covs = paste("V", 1:7, sep = ""),
#'        v_grp = "Group",
#'        cur_grp_level = "current")
#' ps_borrow <- psrwe_borrow(total_borrow = 30, dta_ps)
#' rst       <- psrwe_compl(ps_borrow, v_outcome = "Y_Con")
#' rst
#'
#' @export
#'
psrwe_compl <- function(dta_psbor, v_outcome = "Y",
                      outcome_type = c("continuous", "binary"),
                      ...) {

    ## check
    stopifnot(inherits(dta_psbor,
                       what = get_rwe_class("PSDIST")))

    outcome_type <- match.arg(outcome_type)

    stopifnot(v_outcome %in% colnames(dta_psbor$data))

    ## observed
    rst_obs <- get_observed(dta_psbor$data, v_outcome)

    ## call estimation
    rst <- get_ps_cl_km(dta_psbor, v_outcome = v_outcome,
                        outcome_type = outcome_type,
                        f_stratum = get_cl_stratum, ...)

    ## return
    rst$Observed <- rst_obs
    rst$Method   <- "ps_cl"
    rst$Outcome_type <- outcome_type
    class(rst)   <- get_rwe_class("ANARST")

    rst
}


#' Composite Likelihood Estimation
#'
#' Estimate parameter of interest based composite likelihood for a single PS
#' stratum
#'
#' @inheritParams psrwe_powerp
#'
#' @param dta_cur Vector of outcome from a PS stratum in current study
#' @param dta_ext Vector of outcome from a PS stratum in external data source
#' @param n_borrow Number of subjects to be borrowed
#' @param equal_sd Boolean. whether sd is the same between the current study and
#'     external data source
#'
#' @return Maximum composite likelihood estimator of the mean
#'
#' @examples
#' x <- rnorm(100,  mean = 0, sd = 1)
#' y <- rnorm(1000, mean = 1, sd = 2)
#' rwe_cl(x, y, n_borrow = 20, equal_sd = FALSE)
#'
#' @export
#'
rwe_cl <- function(dta_cur, dta_ext, n_borrow = 0,
                   outcome_type = c("continuous", "binary"),
                   equal_sd = TRUE) {

    f_ll <- function(pars) {
        theta  <- pars[1]
        sig2_1 <- pars[2]
        sig2_0 <- pars[3]

        ll <- - n1 * log(sig2_1) / 2
        ll <- ll - n1 * mean((dta_cur - theta)^2) / 2 / sig2_1
        ll <- ll - n_borrow * log(sig2_0) / 2
        ll <- ll - n_borrow * mean((dta_ext - theta)^2) / 2 / sig2_0

        ll
    }

    f_gradient <- function(pars) {
        theta  <- pars[1]
        sig2_1 <- pars[2]
        sig2_0 <- pars[3]

        g <- numeric(length(pars))

        ## d logl / d theta
        g[1] <- n1 / sig2_1 * (mean(dta_cur) - theta) +
            n_borrow / sig2_0 * (mean(dta_ext) - theta)

        ## d logl / d sig2.1
        g[2] <- - n1 / 2 / sig2_1 +
            n1 * mean((dta_cur - theta)^2) / 2 / sig2_1 / sig2_1
        ## d logl / d sig2.0
        g[3] <- - n_borrow / 2 / sig2_0 +
            n_borrow * mean((dta_cur - theta)^2) / 2 / sig2_0 / sig2_0

        return(g)
    }


    type <- match.arg(outcome_type)
    n1   <- length(dta_cur)

    ## ignore external data
    if (0 == n_borrow) {
        ## placeholder
        dta_ext  <- dta_cur
        equal_sd <- TRUE;
    }

    init_theta <- (n1 / (n1 + n_borrow)) * mean(dta_cur) +
        (n_borrow / (n1 + n_borrow)) * mean(dta_ext)

    if (("continuous" == type & equal_sd) |
        "binary" == type) {
        rst <- init_theta
    } else {
        init_sig2_1 <- mean((dta_cur - init_theta)^2)
        init_sig2_0 <- mean((dta_ext - init_theta)^2)
        rst         <- optim(c(init_theta, init_sig2_1, init_sig2_0),
                             method  = "L-BFGS-B",
                             fn      = f_ll,
                             lower   = c(-Inf, 1e-6, 1e-6),
                             upper   = rep(Inf, 3),
                             control = list(fnscale = -1))$par[1];
    }

    rst
}


#' Get estimation for each stratum
#'
#'
#' @noRd
#'
get_cl_stratum <- function(d1, d0 = NULL, n_borrow = 0, outcome_type, ...) {

    ## treatment or control only
    dta_cur <- d1
    ns1     <- length(dta_cur)
    if (0 == n_borrow | is.null(d0)) {
        theta    <- mean(dta_cur)
        sd_theta <- switch(outcome_type,
                           continuous = sd(dta_cur),
                           binary     = sqrt(theta * (1 - theta) / ns1))

        return(c(theta, sd_theta))
    }

    ## overall ps-cl
    dta_ext <- d0
    ns0     <- length(dta_ext)

    ##  overall estimate
    overall_theta  <- rwe_cl(dta_cur, dta_ext, n_borrow, ...)

    ##jackknife
    jk_theta <- NULL
    for (j in seq_len(ns1)) {
        cur_jk   <- rwe_cl(dta_cur[-j], dta_ext, n_borrow, ...)
        jk_theta <- c(jk_theta, cur_jk)
    }

    if (ns0 > 0) {
        for (j in seq_len(ns0)) {
            cur_jk <- rwe_cl(dta_cur, dta_ext[-j], n_borrow, ...)
            jk_theta <- c(jk_theta, cur_jk)
        }
    }

    ## summary
    sd_theta <- get_jk_sd(overall_theta, jk_theta)
    return(c(overall_theta, sd_theta))
}
