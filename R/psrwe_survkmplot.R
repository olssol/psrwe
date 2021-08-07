#' PS-Integrated Kaplan-Meier Estimation Plot
#'
#' Plot the mean of a survival outcome at all distinctive time points based on
#' PS-integrated Kaplan-Meier approach. Variance is estimated by usual Greenwood
#' method. Applies to the case when there is only one external data source.
#'
#' @inheritParams rwe_ps_powerp
#'
#' @param v_time Column name corresponding to event time
#' @param v_event Column name corresponding to event status
#'
#' @param ... Additional Parameters.
#'
#' @return A data frame with class name \code{RWE_PS_RST_PLOT}. It contains the
#'     composite estimation of the mean for overall and each stratum.
#'
#'
#' @examples
#' \donttest{
#' data(ex_dta)
#' dta_ps <- rwe_ps_est(ex_dta,
#'                      v_covs = paste("V", 1:7, sep = ""),
#'                      v_grp = "Group",
#'                      cur_grp_level = "current")
#' ps_borrow <- rwe_ps_borrow(total_borrow = 30, dta_ps)
#' rst       <- rwe_ps_survkmplot(ps_borrow, v_time = "Y_time",
#'                                v_event = "Status")
#' rst_ci    <- rwe_ps_survkmci(rst)
#' rst_cip   <- rwe_ps_survkmci(rst, conf_type = "plain")
#' plot_survkm(rst_ci)
#' plot_survkm(rst_cip)
#'
#' @export
#'
rwe_ps_survkmplot <- function(dta_psbor,
                              v_time     = "time",
                              v_event    = "event",
                              conf_int = 0.95,
                              conf_type = c("log-log", "plain"),
                              ...) {

    ## check
    stopifnot(inherits(dta_psbor,
                       what = get_rwe_class("PSDIST")))

    stopifnot(all(c(v_event, v_time) %in%
                  colnames(dta_psbor$data)))

    pred_tp <- sort(unique(c(0, dta_psbor$data[[v_time]])))

    ## observed
    rst_obs <- get_km_observed(dta_psbor$data, v_time, v_event, pred_tp)

    ## call estimation
    rst <- get_ps_cl_km(dta_psbor, v_event = v_event, v_time = v_time,
                        f_stratum = get_surv_stratum, pred_tp = pred_tp,
                        f_overall_est = get_overall_est_km,
                        ...)
    ## return
    rst$Observed <- rst_obs
    rst$pred_tp  <- pred_tp
    rst$Method   <- "ps_km"
    class(rst)   <- get_rwe_class("ANARSTPLT")
    rst
}


#' @title S3 method to print RWE_PS_RST_PLOT
#'
#' @method print RWE_PS_RST_PLOT
#'
#' @export
#'
print.RWE_PS_RST_PLOT <- function(x, ...) {
    cat("This is a RWE_PS_RST_PLOT object. Use str() to see details.\n")
    invisible()
}


#' @title S3 method to summary RWE_PS_RST_PLOT
#'
#' @method summary RWE_PS_RST_PLOT
#'
#' @export
#'
summary.RWE_PS_RST_PLOT <- function(object, ...) {
    cat("This is a RWE_PS_RST_PLOT object. Use str() to see details.\n")
    invisible()
}


#' @title Obtain CI for survkm plot.
#'
#' @param dta_pskmp data of RWE_PS_RST_PLOT class
#' @param conf_int level of confidence interval
#' @param conf_type type of confidence interval (log-log as default)
#'
#' @return A data frame with class name \code{RWE_PS_RST_PLOT}. It contains the
#'     composite estimation of the mean for overall and each stratum.
#'     It also contains confidence interval for overall mean.
#'
#' @export
#'
rwe_ps_survkmci <- function(dta_pskmp,
                            conf_int = 0.95,
                            conf_type = c("log-log", "plain")) {

    stopifnot(inherits(dta_pskmp,
                       what = get_rwe_class("ANARSTPLT")))

    ## prepare data
    is_rct <- dta_pskmp$is_rct

    ## get ci
    if (is_rct) {
        ## TBD
    } else {
        S <- dta_pskmp$Control$Overall_Estimate$Mean
        S_se <- dta_pskmp$Control$Overall_Estimate$SD
        rst_km_ci <- get_km_ci(S, S_se,
                               conf_int = conf_int,
                               conf_type = conf_type[1])
        dta_pskmp$Control$Overall_Estimate$CI <- rst_km_ci
    }

    return(dta_pskmp)
}


#' @title Plot KM at all time points
#'
#' @param dta_pskmp data of RWE_PS_RST_PLOT class
#' @param xlab xlab
#' @param ylab ylab
#' @param ylim ylim
#' @param ... Additional Parameters.
#'
#' @return A KM plot.
#'
#' @export
#'
plot_survkm <- function(dta_pskmp,
                        xlab = "Time",
                        ylab = "Survival Probability",
                        ...) {

    stopifnot(inherits(dta_pskmp,
                       what = get_rwe_class("ANARSTPLT")))

    ## prepare data
    is_rct <- dta_pskmp$is_rct
    x      <- dta_pskmp$pred_tp

    ## check arguments
    args <- list(...)
    if ("xlim" %in% names(args)) {
        xlim <- args[['xlim']]
    } else {
        xlim <- range(x)
    }

    ## plot
    if (is_rct) {
        ## TBD
    } else {
        y <- dta_pskmp$Control$Overall_Estimate$Mean
        lower <- dta_pskmp$Control$Overall_Estimate$CI$lower
        upper <- dta_pskmp$Control$Overall_Estimate$CI$upper
        if ("ylim" %in% names(args)) {
            ylim <- args[['ylim']]
        } else {
            ylim <- c(min(lower), max(upper))
        }

        tmp_dta <- data.frame(x = x, y = y, lower = lower, upper = upper)
        rst <- ggplot() +
               geom_step(tmp_dta, mapping = aes(x, y)) +
               geom_step(tmp_dta, mapping = aes(x, lower), linetype = 3) +
               geom_step(tmp_dta, mapping = aes(x, upper), linetype = 3) +
               scale_y_continuous(limits = ylim) +
               scale_x_continuous(limits = xlim) +
               labs(x = xlab, y = ylab)
    }

    return(rst)
}


#' Get KM CI
#'
#' @noRd
#'
get_km_ci <- function(S, S_se,
                      conf_int = 0.95,
                      conf_type = c("log-log", "plain")) {
    z_alphad2 <- qnorm((1 - conf_int) / 2, lower.tail = FALSE)

    if (conf_type[1] == "log-log") {
        log_S <- log(S)
        log_log_S <- log(-log_S)
        # var_log_S <- S_se^2 / S^2
        # se_log_log_S <- sqrt(var_log_S / log_S^2)
        se_log_log_S <- S_se / S / log_S
        A <- cbind(-z_alphad2 * se_log_log_S,
                   z_alphad2 * se_log_log_S)
        ci <- S^exp(A)
    } else if (conf_type[1] == "plain") {
        ci <- cbind(S - z_alphad2 * S_se, S + z_alphad2 * S_se)
    } else {
        stop("conf_type is not implemented.")
    }
    colnames(ci) <- c("lower", "upper")

    rst <- list(conf_int = conf_int,
                conf_type = conf_type,
                lower = ci[, 1],
                upper = ci[, 2])

    return(rst)
}
