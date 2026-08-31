#' Get number of subjects borrowed from each statum
#'
#' Based on PS distances or number of current control subjects, split the total
#' number of subjects to be borrowed from the external data source for each
#' stratum
#'
#' @param dtaps  A class \code{PSRWE_DTA} or \code{PSRWE_DTA_MAT} object.
#' @param total_borrow Total number of subjects to be borrowed
#' @param method Method to split \code{total_borrow} for a class
#'     \code{PSRWE_DTA} object, which can be based on distance (\code{method =
#'     "distance"}) or inverse distance (\code{method = "inverse_distance"}).
#'     Other possible options include \code{"n_current"} and \code{"n_external"}
#'     that use the proportion of stratum sample size based on the current and
#'     external data, respectively.
#'     Ignored for class \code{PSRWE_DTA_MAT} object.
#' @param .drop_arg_fml internal use to drop arguments and call, this is
#'     only used in cjk.
#' @param ... Additional parameters for \code{\link{summary.PSRWE_DTA}}.
#'
#' @return A class \code{PSRWE_BORR} list. It appends the following items to
#'     the \code{dtaps}:
#'
#' \describe{
#'
#'     \item{Proportion}{Proportion splitting the number of total borrow among
#'     strata.}
#'
#'     \item{N_Borrow}{The number of subjects to be borrowed in each stratum.}
#'
#'     \item{Alpha}{Weight parameter value in each stratum.}
#' }
#'
#' @examples
#' data(ex_dta)
#' dta_ps <- psrwe_est(ex_dta,
#'                      v_covs = paste("V", 1:7, sep = ""),
#'                      v_grp = "Group",
#'                      cur_grp_level = "current")
#'
#' ps_borrow <- psrwe_borrow(total_borrow = 20, dta_ps)
#' ps_borrow
#'
#' ## Use different similarity metric
#' ps_borrow_omkss <- psrwe_borrow(total_borrow = 20, dta_ps,
#'                                  metric = "omkss")
#' ps_borrow_omkss
#'
#' @export
#'
#'
psrwe_borrow <- function(dtaps, total_borrow,
                          method = c("distance",
                                     "inverse_distance",
                                     "n_current",
                                     "n_external"),
                          .drop_arg_fml = FALSE,
                          ...) {

    ## save arguments and call first (overwrite the drop if in cjk)
    if (.drop_arg_fml) {
        call_arg <- NA
        call_fml <- NA
    } else {
        call_arg <- c(as.list(environment()), list(...))
        call_arg[["dtaps"]] <- NA
        call_fml <- as.character(match.call()[[1]])
    }

    ## prepare
    is_ps       <- inherits(dtaps, what = get_rwe_class("DWITHPS"))
    is_ps_match <- inherits(dtaps, what = get_rwe_class("DPSMATCH"))
    stopifnot(is_ps | is_ps_match)

    method   <- match.arg(method)

    rst_sum <- summary(dtaps, ...)
    ns0     <- rst_sum$Summary$N_RWD
    rs      <- rst_sum$Summary$Distance

    if (dtaps$is_rct) {
        ns1 <- rst_sum$Summary$N_Cur_CTL
    } else {
        ns1 <- rst_sum$Summary$N_Current
    }

    if (is_ps_match) {
        method <- "n_current"
        ## TODO: The distance may base on matched samples
    }

    borrow  <- get_aborrow(total_borrow,
                           ns0, ns1, rs,
                           m_lambda = method, ...)

    ## return
    dtaps$Total_borrow  <- total_borrow
    dtaps$Borrow        <- cbind(rst_sum$Summary, borrow)
    dtaps$Borrow_method <- method
    dtaps$Call_arg$psrwe_borrow <- call_arg
    dtaps$Call_fml$psrwe_borrow <- call_fml

    class(dtaps) <- get_rwe_class("PSDIST")
    return(dtaps)
}

#' @title Print borrow information
#'
#' @description Print summary information of borrowing
#'
#' @param x A list of class \code{PSRWE_BOR} that is generated using
#'   the \code{\link{psrwe_borrow}} function.
#'
#' @param ... Additional parameters
#'
#' @seealso  \code{\link{psrwe_borrow}}
#'
#' @method print PSRWE_BOR
#'
#' @return A list from \code{x$Borrow}
#'
#' @export
#'
print.PSRWE_BOR <- function(x, ...) {
    ss <- paste("A total of ",
                x$Total_borrow,
                " subjects will be borrowed from the RWD. ",
                "The number ",
                x$Total_borrow,
                " is split proportional to ",
                get_rwe_class(x$Borrow_method),
                " in each stratum. The following table",
                " summarizes the number of subjects to be borrowed",
                " and the weight parameter in each stratum:",
                sep = "")
    cat_paste(ss)
    cat("\n")

    print(x$Borrow)
}
