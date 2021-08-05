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
#' @return A data frame with class name \code{RWE_PS_RST}. It contains the
#'     composite estimation of the mean for overall and each stratum.
#'
#'
#' @examples
#' \donttest{
#' data(ex_dta)
#' dta_ps <- rwe_ps_est(ex_dta,
#'        v_covs = paste("V", 1:7, sep = ""),
#'        v_grp = "Group",
#'        cur_grp_level = "current")
#' ps_borrow <- rwe_ps_borrow(total_borrow = 30, dta_ps)
#' rst_km_allt <- rwe_ps_survkmplot(ps_borrow, v_time = "Y_time",
#'                                  v_event = "Status")
#' plot_survkm(rst_km_allt,, ylim = c(0.4, 1))
#' plot
#'
#' @export
#'
rwe_ps_survkmplot <- function(dta_psbor,
                              v_time     = "time",
                              v_event    = "event",
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
#' @noRd
#'
summary.RWE_PS_RST_PLOT <- function(object, ...) {
  cat("This is a RWE_PS_RST_PLOT object. Use str() to see details.\n")
  invisible()
}


#' @title Plot KM at all time points
#'
#' @param dta_pskmp data
#' @param conf.int confidence interval
#' @param conf.type confidence type (log as default)
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
                        conf.int = 0.95,
                        conf.type = c("log"),
                        xlab = "Time",
                        ylab = "Survival Probability",
                        ylim = c(0, 1),
                        ...) {

    stopifnot(inherits(dta_pskmp,
                       what = get_rwe_class("ANARSTPLT")))

    ## prepare data
    is_rct <- dta_pskmp$is_rct
    x      <- dta_pskmp$pred_tp

    ## check arguments
    args <- list(...)
    if ("xlim" %in% names(args)) {
      xlim <- args['xlim']
    } else {
      xlim <- range(x)
    }

    if (is_rct) {
      ## TBD
    } else {
      y <- dta_pskmp$Control$Overall$Mean
      tmp.dta <- data.frame(x = x, y = y) 
      rst <- ggplot() +
             geom_step(data = tmp.dta, mapping = aes(x = x, y = y)) +
             scale_y_continuous(limits = ylim) +
             scale_x_continuous(limits = xlim) +
             labs(x = xlab, y = ylab)
    }

    return(rst)
}
