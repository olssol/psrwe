#' @title Estimate propensity scores
#'
#' @description Estimate propensity scores using logistic regression or random
#'   forest model.
#'
#' @param data Data frame with group assignment and covariates.
#' @param ps_fml Propensity score (PS) formula. If \code{NULL}, all covariates
#'     will be included in the PS model in a linear form.
#' @param ps_method Method for PS estimation. \code{logistic}: By logistic
#'     regression; \code{radomforest}: By randomForest model.
#' @param v_covs Column names corresponding to covariates.
#' @param v_grp Column name corresponding to group assignment.
#' @param cur_grp_level Group level for the current study. Default is
#'     \code{cur_grp_level = 1}. Ignored for single-arm studies.
#' @param v_arm Column name corresponding to arm assignment.
#' @param ctl_arm_level Arm level for the control arm. Ignored for single-arm
#'     studies.
#' @param stra_ctl_only Create strata by control arm patients only. Default
#'     \code{TRUE}. Ignored by single-arm studies. For randomized studies, when
#'     \code{stra_ctl_only} is \code{FALSE}, strata are created based on the PS
#'     scores of the entire current study patients.
#' @param nstrata Number of PS strata to be created.
#' @param ps_method Method to calculate propensity scores. Can be set to
#'     \code{logistic} for logistic regression or \code{randomforest} for a
#'     random forest approach.
#' @param trim_ab Trim external subjects who are above or below the
#'     range of current study. Default \code{both} trims both above and below.
#'     Other options include \code{above} for above only, \code{below} for
#'     below only, and \code{none} for no trimming.
#' @param .drop_arg_fml internal use to drop arguments and call, this is
#'     only used in cjk.
#' @param ... Additional parameters for calculating the propensity score to be
#'     used in \code{randomForest} or \code{glm} .
#'
#' @return A list of class \code{PSRWE_DAT} with items:
#'
#' \describe{
#'   \item{data}{Original data with column \code{_ps_} for estimated PS scores
#'               and \code{_strata_} for PS stratum added.}
#'   \item{ps_fml}{PS formula for estimated PS scores.}
#'   \item{is_rct}{Whether the current study is a randomized study.}
#'   \item{nstrata}{Number of strata.}
#' }
#'
#' @examples
#' data(ex_dta)
#' psrwe_est(ex_dta,
#'        v_covs = paste("V", 1:7, sep = ""),
#'        v_grp = "Group",
#'        cur_grp_level = "current")
#'
#' @export
#'
psrwe_est <- function(data,
                       ps_fml        = NULL,
                       ps_method     = c("logistic", "randomforest"),
                       v_covs        = "V1",
                       v_grp         = "Group",
                       cur_grp_level = 1,
                       v_arm         = NULL,
                       ctl_arm_level = NULL,
                       stra_ctl_only = TRUE,
                       nstrata       = 5,
                       trim_ab       = c("both", "above", "below", "none"),
		       .drop_arg_fml = FALSE,
                       ...) {

    ## save arguments and call first (overwrite the drop if in cjk)
    if (.drop_arg_fml) {
        call_arg <- NA
        call_fml <- NA
    } else {
        call_arg <- c(as.list(environment()), list(...))
        call_arg[["data"]] <- NA
        call_fml <- as.character(match.call()[[1]])
    }

    ## prepare
    if (!identical("data.frame", class(data))) {
        warning("data should be a data.frame object")
    }

    data      <- as.data.frame(data)
    ps_method <- match.arg(ps_method)
    trim_ab   <- match.arg(trim_ab)

    ## Generate formula
    if (is.null(ps_fml)) {
        ps_fml <- as.formula(paste(v_grp, "~",
                                   paste(v_covs, collapse = "+"),
                                   sep = ""))
    }

    dnames <- colnames(data)
    v_fml  <- all.vars(ps_fml)
    if (!(all(v_fml %in% dnames))) {
        stop("Group or covariate was not found in data")
    }

    ## assign v_grp if ps_fml was used
    v_grp <- v_fml[1]

    ## Arm
    if (!is.null(v_arm)) {
        if (!(v_arm %in% dnames))
            stop("Arm was not found in data")

        is_rct <- TRUE
    } else {
        is_rct <- FALSE
    }

    # Current study index will be kept in the results
    d1_inx   <- cur_grp_level == data[[v_grp]]
    keep_inx <- which(d1_inx)
    if (0 == length(keep_inx))
        stop("No current study subjects found in data")

    # For two-arm studies only
    if (!is.null(v_arm) &
        stra_ctl_only) {
        d1_inx <- d1_inx & ctl_arm_level == data[[v_arm]]
    }

    # Get propensity scores
    all_ps <- get_ps(data, ps_fml = ps_fml, ps_method = ps_method, ...)

    # Add groups (0/1), arms (0/1), ps to data
    data[["_ps_"]]  <- all_ps

    grp <- rep(1, nrow(data))
    grp[which(data[[v_grp]] != cur_grp_level)] <- 0
    data[["_grp_"]] <- grp

    arm <- rep(0, nrow(data))
    if (is_rct) {
        arm[which(data[[v_arm]] != ctl_arm_level)] <- 1
    }
    data[["_arm_"]]  <- arm

    # Stratification
    if (nstrata > 0) {
        d1_ps   <- all_ps[which(d1_inx)]
        strata  <- rwe_cut(d1_ps, all_ps,
                           breaks   = nstrata,
                           keep_inx = keep_inx,
                           trim_ab  = trim_ab)

        data[["_strata_"]] <- factor(strata,
                                     levels = c(1:nstrata),
                                     labels = paste("Stratum", 1:nstrata))
    }

    # Add new ID
    data <- data %>%
        arrange(`_grp_`, `_arm_`, `_strata_`) %>%
        mutate(`_id_` = 1:n())

    # Outputx
    rst <- list(data      = data,
                ps_fml    = ps_fml,
                ps_method = ps_method,
                is_rct    = is_rct,
                nstrata   = nstrata,
                trim_ab   = trim_ab,
                Call_arg  = list(psrwe_est = call_arg),
                Call_fml  = list(psrwe_est = call_fml))

    class(rst) <- get_rwe_class("DWITHPS")
    return(rst)
}

#' @title Summarize PS estimation and stratification results
#'
#' @description Get the number of subjects and the distances of PS
#'   distributions for each PS stratum.
#'
#' @inheritParams get_distance
#'
#' @param object A list of class \code{PSRWE_DAT} that is generated using
#'   the \code{\link{psrwe_est}} function.
#' @param min_n0 threshold for the number of external subjects, below which the
#'   external data in the current stratum will be ignored by setting the PS
#'   distance to 0. Default value 10.
#' @param ... Additional parameters.
#'
#' @return A list with columns:
#'
#' \describe{
#'     \item{Summary}{A data frame with Stratum, number of subjects in RWD,
#'     current study, number of subjects in control and treatment arms for RCT
#'     studies, and distance in PS distributions.}
#'
#'     \item{Overall}{A data frame with the overall number of not-trimmed
#'     subjects in RWD, number of patients in the current study, number of
#'     subjects in control and treatment arms for RCT studies, and distance
#'     in PS distributions.}
#'
#'     \item{N}{Vector of total number of total RWD patients, number of trimmed
#'     RWD patients, and total number of current study patients.}
#'
#'     \item{ps_fml}{PS model.}
#'     \item{Distance_metric}{Metric used for calculating the distance.}
#' }
#'
#' @method summary PSRWE_DTA
#'
#' @examples
#' data(ex_dta)
#' dta_ps <- psrwe_est(ex_dta,
#'                      v_covs = paste("V", 1:7, sep = ""),
#'                      v_grp = "Group",
#'                      cur_grp_level = "current")
#' dta_ps
#'
#' ## With different similarity metric
#' print(dta_ps, metric = "omkss")
#' dta_ps_sum <- summary(dta_ps, metric = "omkss")
#'
#' @export
#'
summary.PSRWE_DTA <- function(object,
                                metric = c("ovl", "ksd", "std", "abd",
                                           "ley", "mhb", "omkss"),
                                min_n0 = 10, ...) {

    f_narm <- function(inx) {
      if (is_rct) {
          n0  <- length(which(0 == dataps[inx, "_arm_"]))
          n1  <- length(which(1 == dataps[inx, "_arm_"]))
          rst <- c(length(inx), n0, n1)
      } else {
          rst <- length(inx)
      }
      rst
    }

    metric  <- match.arg(metric)
    dataps  <- object$data
    nstrata <- object$nstrata
    is_rct  <- object$is_rct
    strata  <- levels(dataps[["_strata_"]])
    trim_ab <- object$trim_ab

    if (is_rct) {
        col_n  <- c("N_RWD",     "N_RWD_CTL", "N_RWD_TRT",
                    "N_Current", "N_Cur_CTL", "N_Cur_TRT",
                    "Distance")
    } else {
        col_n  <- c("N_RWD", "N_Current", "Distance")
    }

    n_trim    <- length(which(is.na(dataps[["_strata_"]])))
    n_rwd     <- length(which(0 == dataps[["_grp_"]]))
    n_current <- nrow(dataps) - n_rwd

    ps_range_current <- range(dataps[["_ps_"]][dataps[["_grp_"]] == 1])
    n_rwd_below <- sum(dataps[["_ps_"]][dataps[["_grp_"]] == 0] <
                       ps_range_current[1])
    n_rwd_above <- sum(dataps[["_ps_"]][dataps[["_grp_"]] == 0] >
                       ps_range_current[2])

    rst <- NULL
    for (i in strata) {
        inx_ps0 <- i == dataps[["_strata_"]] & 0 == dataps[["_grp_"]]
        inx_ps1 <- i == dataps[["_strata_"]] & 1 == dataps[["_grp_"]]
        n0_01   <- f_narm(which(inx_ps0))
        n1_01   <- f_narm(which(inx_ps1))

        if (is_rct) {
            inx_ps0 <- inx_ps0 & 0 == dataps[["_arm_"]]
            inx_ps1 <- inx_ps1 & 0 == dataps[["_arm_"]]
        }

        ps0 <- dataps[which(inx_ps0), "_ps_"]
        ps1 <- dataps[which(inx_ps1), "_ps_"]

        if (0 == length(ps0) | 0 == length(ps1)) {
            warning(paste("No samples in", i))
        }

        if (any(is.na(c(ps0, ps1)))) {
            warning(paste("NA found in propensity scores in ", i))
        }

        if (length(ps0) < min_n0) {
            warning(paste("Not enough data in the external data in ",
                          i, ". External data ignored.",
                          sep = ""))

            cur_dist <- 0
        } else {
            cur_dist <- get_distance(ps0, ps1, metric[1])
        }

        rst <- rbind(rst, c(n0_01, n1_01, cur_dist))
    }
    colnames(rst) <- col_n

    ## overall
    inx_tot_ps0 <- which(0 == dataps[["_grp_"]] &
                         !is.na(dataps[["_strata_"]]))
    inx_tot_ps1 <- which(1 == dataps[["_grp_"]])
    n0_tot_01   <- f_narm(inx_tot_ps0)
    n1_tot_01   <- f_narm(inx_tot_ps1)

    ps0         <- dataps[inx_tot_ps0, "_ps_"]
    ps1         <- dataps[inx_tot_ps1, "_ps_"]
    all_dist    <- get_distance(ps0, ps1, metric[1])

    rst_overall           <- rbind(c(n0_tot_01, n1_tot_01, all_dist))
    colnames(rst_overall) <- col_n

    # Return
    rst <- list(Summary         = cbind(data.frame(Stratum = strata),
                                        data.frame(rst)),
                Overall         = cbind(data.frame(Stratum = "Overall"),
                                        rst_overall),
                N               = c(RWD     = n_rwd,
                                    Current = n_current,
                                    Trimmed = n_trim,
                                    RWD_below_current = n_rwd_below,
                                    RWD_above_current = n_rwd_above,
                                    Trim_ab = trim_ab),
                ps_fml          = object$ps_fml,
                Distance_metric = metric[1])

    invisible(rst)
}

#' @title Print PS estimation results
#'
#' @description Print summary information of PS estimation results
#'
#' @param x A list of class \code{PSRWE_DTA} that is generated using
#'   the \code{\link{psrwe_est}} function.
#' @param ... Parameters for \code{summery.PSRWE_DTA}
#'
#' @seealso  \code{\link{summary.PSRWE_DTA}}
#'
#'
#' @method print PSRWE_DTA
#'
#' @return A list from \code{summary(x)} with additional information
#'
#' @export
#'
print.PSRWE_DTA <- function(x, ...) {
    rst_sum <- summary(x, ...)

    ## overall summary
    cat_ps_dta(x, rst_sum)

    ## summary table
    cat("\n")
    ss <- paste("The following table summarizes the number of subjects ",
                "in each stratum, and the ",
                "distance in PS distributions ", "calculated by ",
                get_rwe_class(rst_sum$Distance_metric), ":", sep = "")

    cat_paste(ss)
    cat("\n")

    ## number of patients
    tbl_summary           <- rst_sum$Summary
    tbl_summary$N_RWD_TRT <- NULL
    print(tbl_summary)
}


#' @title Plot PS distributions
#'
#' @description S3 method for visualizing PS adjustment
#'
#' @param x Class \code{RWE_DWITHPS} created by \code{psrwe_*} functions
#' @param plot_type Types of plots.
#' \describe{
#'   \item{ps}{PS density plot}
#'   \item{balance}{Covariate balance plot}
#'   \item{diff}{Standardized mean differences, metric = std or astd}
#' }
#' @param ... Additional parameter for the plot
#'
#' @method plot PSRWE_DTA
#'
#' @return A plot of class in ggplot2
#'
#' @export
#'
plot.PSRWE_DTA <- function(x, plot_type = c("ps", "balance", "diff"), ...) {
    type <- match.arg(plot_type)
    switch(type,
           ps      = plot_ps(x, ...),
           balance = plot_balance(x, ...),
           diff    = plot_astd(x, ...))
}


#' @title Create strata
#'
#' @description
#' Cut a sequence of numbers into bins.
#'
#' The cut points are chosen such that there will be equal numbers in each bin
#' for \code{x}. By default, values of \code{y} that are outside the range of
#' \code{x} will be excluded from the bins, unless they are in the
#' \code{keep_inx}.
#'
#' @param x Vector of values based on which cut points will be determined
#' @param y Vector of values to be cut, default to be the same as \code{x}
#' @param breaks Number of cut points
#' @param keep_inx Indices of y that will be categorized as 1 or the largest bin
#'     even if their values are out of range of x, i.e. the y's that will not be
#'     trimmed
#' @param trim_ab Trim external subjects who are above or below the
#'     range of current study. Default \code{both} trims both above and below.
#'     Other options include \code{above} for above only, \code{below} for
#'     below only, and \code{none} for no trimming.
#'
#' @return A vector of stratum assignment for \code{y}. The y's that are outside
#'     the range of \code{x} and not in \code{keep_inx} are assigned \code{NA}
#'     in the result.
#'
#' @examples
#'
#' x <- rnorm(100,  mean = 0, sd = 1)
#' y <- rnorm(1000, mean = 1, sd = 2)
#' rwe_cut(x, y, breaks = 5)
#'
#' @export
#'
#'
rwe_cut <- function(x, y = x, breaks = 5, keep_inx = NULL,
                    trim_ab = c("both", "above", "below", "none")) {
    cuts    <- quantile(x, seq(0, 1, length = breaks + 1))
    cuts[1] <- cuts[1] - 1e-100

    ## Overwrite cuts[1] and cuts[breaks]
    trim_ab <- match.arg(trim_ab)
    switch(trim_ab,
           none = {
             cuts[1] <- -Inf
             cuts[breaks + 1] <- Inf
           },
           above = {
             cuts[1] <- -Inf
           },
           below = {
             cuts[breaks + 1] <- Inf
           })

    rst     <- rep(NA, length(y))
    for (i in 2:length(cuts)) {
        inx      <- which(y > cuts[i - 1] & y <= cuts[i])
        rst[inx] <- i - 1
    }

    if (!is.null(keep_inx)) {
        inx <- which(y[keep_inx] <= cuts[1])
        if (0 < length(inx)) {
            rst[keep_inx[inx]] <- 1
        }

        inx <- which(y[keep_inx] > cuts[length(cuts)])
        if (0 < length(inx)) {
            rst[keep_inx[inx]] <- length(cuts) - 1
        }
    }

    rst
}
