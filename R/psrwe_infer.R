#' @title Inference for the PS-Integrated Estimation
#'
#' @description
#' Inference for the PS-integrated approach.
#'
#' @param dta_psrst A returned object with class \code{PSRWE_EST}
#' @param alternative A character string for the alternative hypothesis that
#'        must be one of \code{"less"} (default), \code{"greater"}, or
#'        \code{"two_sided"} (for log-rank and RMST only)
#' @param mu A number indicating the true value of the parameter of interest
#'        (or the difference in means for two arms),
#'        \code{mu = 0} when the test is log-rank or RMST
#' @param method_pval A method name for p-value (default wald),
#'        no effect on Bayesian methods, and
#'        \code{method = "score"} only is for binary outcomes in
#'        single arm studes (i.e., comparing with a PG set by \code{mu})
#' @param ... Other options
#'
#' @return A list with class name \code{PSRWE_EST}.
#'
#' @examples
#' data(ex_dta)
#' dta_ps <- psrwe_est(ex_dta,
#'        v_covs = paste("V", 1:7, sep = ""),
#'        v_grp = "Group",
#'        cur_grp_level = "current")
#' ps_borrow <- psrwe_borrow(total_borrow = 30, dta_ps)
#' ps_rst <- psrwe_compl(ps_borrow, v_outcome = "Y_Con")
#' rst <- psrwe_infer(ps_rst)
#' rst
#'
#' @export
#'
psrwe_infer <- function(dta_psrst,
                        alternative = c("less", "greater", "two_sided"),
                        mu = 0,
                        method_pval = c("wald", "score", "score_weighted"),
                        ...) {

    ## check
    stopifnot(inherits(dta_psrst,
                       what = get_rwe_class("ANARST")))

    stopifnot(dta_psrst$Method %in% get_rwe_class("ANAMETHOD"))

    alternative <- match.arg(alternative)
    method_pval <- match.arg(method_pval)

    ## get pval by method
    if (dta_psrst$Method == "ps_pp") {
        if (alternative == "two_sided") {
             stop("two_sided is not implemented for ps_pp")
        }
        rst_psinfer <- get_psinfer_bayesian(dta_psrst,
                                            alternative,
                                            mu)
    } else {
        if (dta_psrst$Method %in% c("ps_lrk", "ps_rmst")) {
            alternative <- "two_sided"
            mu <- 0
        }
        rst_psinfer <- get_psinfer_freq(dta_psrst,
                                        alternative,
                                        mu,
                                        method_pval)
    }

    ## return
    rst <- dta_psrst
    rst$INFER <- rst_psinfer
    return(rst)
}




#' @title Bayesian inference
#'
#' @noRd
get_psinfer_bayesian <- function(dta_psrst,
                                 alternative,
                                 mu) {

    ## prepare data
    is_rct <- dta_psrst$is_rct 
    if (is_rct) {
        type <- "Effect"
    } else {
        type <- "Control"
    }

    ## prepare for the return object
    rst_psinfer <- list(Control = NULL,
                        Effect = NULL,
                        Method_infer = "posterior probability",
                        Alternative = alternative,
                        Mu = mu,
                        Method_pval = NA)

    ## by study type
    rst_psinfer[[type]]$Stratum_InferProb <-
        get_bpostp(dta_psrst[[type]]$Stratum_Samples,
                   alternative = alternative,
                   mu = mu)
    rst_psinfer[[type]]$Overall_InferProb <-
        get_bpostp(dta_psrst[[type]]$Overall_Samples,
                   alternative = alternative,
                   mu = mu)

    return(rst_psinfer)
}


#' @title Bayesian posterior probability
#'
#' @noRd
get_bpostp <- function(x,
                       alternative,
                       mu) {

    if (is.vector(x)) {
        post_prob <- switch(alternative,
                            less = {
                              mean(x < mu)
                            },
                            greater = {
                              mean(x > mu)
                            })
    } else {
        post_prob <- switch(alternative,
                            less = {
                              rowMeans(x < mu)
                            },
                            greater = {
                              rowMeans(x > mu)
                            })
    }

    rst <- data.frame(Infer_prob = post_prob)
    return(rst)
}




#' @title Frequentist p-value
#'
#' @noRd
get_psinfer_freq <- function(dta_psrst,
                             alternative,
                             mu,
                             method_pval) {

    ## prepare data
    is_rct <- dta_psrst$is_rct 
    if (is_rct) {
        type <- "Effect"
    } else {
        type <- "Control"
    }
    outcome_type <- dta_psrst$Outcome_type

    ## prepare for the return object
    rst_psinfer <- list(Control = NULL,
                        Effect = NULL,
                        Method_infer = "p_value",
                        Alternative = alternative,
                        Mu = mu,
                        Method_pval = method_pval)

    ## by method_pval, study type, and outcome type
    if (method_pval %in% c("score", "score_weighted")) {
      if (is_rct || outcome_type != "binary") {
        stop("socre pval is only for binary outcomes and single arm study.")
      } else {
        N_borrow <- dta_psrst$Borrow$N_Borrow
        N_current <- dta_psrst$Borrow$N_Current
        N_nominal <- N_current + N_borrow

        rst_psinfer[[type]]$Stratum_InferProb <-
            get_fpval_binary_score(dta_psrst[[type]]$Stratum_Estimate,
                                   alternative = alternative,
                                   mu = mu,
                                   n_nominal = N_nominal)

        if (method_pval == "score") {
          rst_psinfer[[type]]$Overall_InferProb <-
              get_fpval_binary_score(dta_psrst[[type]]$Overall_Estimate,
                                     alternative = alternative,
                                     mu = mu,
                                     n_nominal = sum(N_nominal))
        } else {
          rst_psinfer[[type]]$Overall_InferProb <-
              get_fpval_binary_score_weighted(
                  dta_psrst[[type]]$Overall_Estimate,
                  alternative = alternative,
                  mu = mu,
                  n_nominal = N_nominal,
                  weights = N_current)
        }
      }
    } else {
      rst_psinfer[[type]]$Stratum_InferProb <-
          get_fpval(dta_psrst[[type]]$Stratum_Estimate,
                    alternative = alternative,
                    mu = mu)
      rst_psinfer[[type]]$Overall_InferProb <-
          get_fpval(dta_psrst[[type]]$Overall_Estimate,
                    alternative = alternative,
                    mu = mu)
    }

    return(rst_psinfer)
}


#' @title Frequentist p-value (normal approximated)
#'
#' @noRd
get_fpval <- function(x,
                      alternative,
                      mu) {

    tstat <- (x$Mean - mu) / x$StdErr
    p_value <- switch(alternative,
                      less = {
                        pnorm(tstat)
                      },
                      greater = {
                        pnorm(tstat, lower.tail = FALSE)
                      },
                      two_sided = {
                        pnorm(abs(tstat), lower.tail = FALSE) * 2
                      })

    rst <- data.frame(Infer_prob = p_value)
    return(rst)
}


#' @title Frequentist p-value (for binary outcome and score method)
#'
#' @noRd
get_fpval_binary_score <- function(x,
                                   alternative,
                                   mu,
                                   n_nominal) {

    ## mu is the PG0
    stderr <- sqrt(mu * (1 - mu) / n_nominal)
    tstat <- (x$Mean - mu) / stderr
    p_value <- switch(alternative,
                      less = {
                        pnorm(tstat)
                      },
                      greater = {
                        pnorm(tstat, lower.tail = FALSE)
                      },
                      two_sided = {
                        pnorm(abs(tstat), lower.tail = FALSE) * 2
                      })

    rst <- data.frame(Infer_prob = p_value)
    return(rst)
}


#' @title Frequentist p-value (for binary outcome and score_weighted method)
#'
#' @noRd
get_fpval_binary_score_weighted <- function(x,
                                            alternative,
                                            mu,
                                            n_nominal,
                                            weights) {

    ## mu is the PG0
    ws <- weights / sum(weights)
    stderr <- sqrt((mu * (1 - mu) / n_nominal) %*% ws^2)
    tstat <- (x$Mean - mu) / stderr
    p_value <- switch(alternative,
                      less = {
                        pnorm(tstat)
                      },
                      greater = {
                        pnorm(tstat, lower.tail = FALSE)
                      },
                      two_sided = {
                        pnorm(abs(tstat), lower.tail = FALSE) * 2
                      })

    rst <- data.frame(Infer_prob = p_value)
    return(rst)
}

