#' @title PS matching
#'
#' @description Match patients in the external data source with patients in
#'     the current study based on PS using the nearest neighbor method.
#'
#' @param dta_ps A list of class \code{PSRWE_DAT} that is generated using the
#'     \code{\link{psrwe_est}} function.
#' @param ratio Matching ratio (RWD : Current). Default is 3
#'     (i.e., 3:1 matching).
#' @param strata_covs Covariates used for stratification in matching.
#' @param caliper PS matching caliper width. Default 1. This specifies a
#'     width (Euclidean distance) on the probability scale.
#' @param seed Random seed.
#' @param method matching algorithm for PS matching.
#' @param .drop_arg_fml internal use to drop arguments and call, this is
#'     only used in cjk.
#' @param ... Additional parameters for matching
#'
#' @return A list of class \code{PSRWE_DTA_MAT} with items:
#'
#' \describe{
#'
#' \item{data}{Original data with column \code{_ps_} for estimated PS scores,
#'   \code{match_id} for matched current study subject ID, and \code{_strata_}
#'   for PS stratum added.}
#' \item{ps_fml}{PS formula for estimated PS scores.}
#' \item{nstrata}{Number of strata.}
#' \item{ratio}{Matching ratio.}
#'
#' }
#'
#' @examples
#' data(ex_dta)
#' dta_ps <- psrwe_est(ex_dta,
#'                      v_covs = paste("V", 1:7, sep = ""),
#'                      v_grp = "Group",
#'                      cur_grp_level = "current")
#'
#' ## With neighest neighbor
#' .remark_nnmatch <- function() {
#'   dta_ps_mat <- psrwe_match(dta_ps, ratio = 2, strata_covs = "V1",
#'                             seed = 123)
#'   dta_ps_mat
#' }
#' ## Unmark below to run
#' # .remark_nnmatch()
#'
#' ## With optmatch
#' .remark_optmatch <- function() {
#'     warning("The optmatch package may restrict use (academic license).")
#'     dta_ps_opt <- psrwe_match(dta_ps, ratio = 2,
#'                                strata_covs = "V1",
#'                                method = "optm", seed = 123)
#'     dta_ps_opt
#' }
#' ## Unmark below to run if optmatch is available.
#' # .remark_optmatch()
#'
#' @export
#'
psrwe_match <- function(dta_ps, ratio = 3, strata_covs  = NULL,
                         caliper = 1, seed = NULL,
                         method = c("nnwor", "optm"),
			 .drop_arg_fml = FALSE,
                         ...) {

    ## save arguments and call first (overwrite the drop if in cjk)
    if (.drop_arg_fml) {
        call_arg <- NA
        call_fml <- NA
    } else {
        call_arg <- c(as.list(environment()), list(...))
        call_arg[["dta_ps"]] <- NA
        call_fml <- as.character(match.call()[[1]])
    }

    ## prepare
    stopifnot(get_rwe_class("DWITHPS") %in% class(dta_ps))

    mat_method <- match.arg(method)

    ## save the seed from global if any then set random seed
    # old_seed <- NULL
    if (!is.null(seed)) {
        # if (exists(".Random.seed", envir = .GlobalEnv)) {
        #     old_seed <- get(".Random.seed", envir = .GlobalEnv)
        # }
        set.seed(seed)
    }

    ## check stratification factors
    data   <- dta_ps$data
    dnames <- colnames(data)
    if (!(all(strata_covs %in% dnames))) {
        stop("Stratification covariate was not found in data")
    }

    ## get strata
    dat_strata <- get_strata(data, strata_covs)
    nstrata    <- dat_strata$nstrata
    if (nstrata > 5)
        warning("Number of strata is more than 5. Data may be over-stratified.")

    ## prepare data
    data <- dat_strata$data

    ## start matching
    data <- switch(mat_method,
                   nnwor = get_match_nnwor(data, ratio, caliper, ...),
                   optm  = get_match_optm(data, ratio, caliper, ...)
                   )

    ## reset the orignal seed back to the global or
    ## remove the one set within this session earlier.
    # if (!is.null(seed)) {
    #     if (!is.null(old_seed)) {
    #         invisible(assign(".Random.seed", old_seed, envir = .GlobalEnv))
    #     } else {
    #         invisible(rm(list = c(".Random.seed"), envir = .GlobalEnv))
    #     }
    # }

    ## result
    rst             <- dta_ps
    rst$data        <- as.data.frame(data)
    rst$nstrata     <- nstrata
    rst$ratio       <- ratio
    rst$caliper     <- caliper
    rst$strata_covs <- strata_covs
    rst$mat_method  <- mat_method
    rst$Call_arg$psrwe_match <- call_arg
    rst$Call_fml$psrwe_match <- call_fml
    class(rst)      <- get_rwe_class("DPSMATCH")
    return(rst)
}


#' @title Summarize PS estimation and matching results
#'
#' @description Get number of subjects for each PS stratum.
#'
#' @param object A list of class \code{PSRWE_DTA_MAT} that is generated using
#'     the \code{\link{psrwe_match}} function.
#'
#' @param ... Additional parameters.
#'
#' @return A list with columns:
#'
#' \describe{
#'
#'     \item{Summary}{A data frame with Stratum (defined by covariates), number
#'     of subjects in RWD, current study, number of subjects in control and
#'     treatment arms for RCT studies.}
#'
#'     \item{Overall}{A data frame with the overall number of not-trimmed
#'     subjects
#'     in RWD, number of patients in the current study, number of subjects in
#'     control and treatment arms for RCT studies.}
#'
#'     \item{N}{Vector of total number of total RWD patients, number of trimmed
#'     RWD patients, total number of current study patients, number of current
#'     control patients with less than \code{ratio} matched RWD subjects.}
#'
#'     \item{ps_fml}{PS model.}
#'
#'     \item{N_Match}{Number of current control subjects matched with
#'     \code{ratio}, 0 and other number of RWD subjects.}
#'
#'     \item{ratio}{Matching ratio.}
#' }
#'
#' @method summary PSRWE_DTA_MAT
#'
#'
#' @export
#'
summary.PSRWE_DTA_MAT <- function(object, ...) {
    rst_sum <- summary.PSRWE_DTA(object, ...)

    ## adjust rst_sum
    rst_sum$Distance_metric  <- NULL
    rst_sum$Summary$Distance <- NULL
    rst_sum$Overall$Distance <- NULL
    rst_sum$ratio            <- object$ratio

    ## TODO: The distance may base on matched samples
    # if ("metric" %in% names(list(...))) {
        # rst_sum$Summary$Distance <- update_distance(object, ...)
        ## or
        # rst_sum <- update_distance(object, ...)
    # }

    ## check matching ratio
    match_n   <- object$data       %>%
        dplyr::filter(1 == `_grp_` &
                      0 == `_arm_`) %>%
        select(`_matchn_`)

    n_1     <- length(which(object$ratio == match_n))
    n_0     <- length(which(0 == match_n))
    n_other <- nrow(match_n) - n_1 - n_0
    rst_sum$N_Match <- c(Ratio = n_1, Zero = n_0, Other = n_other)

    invisible(rst_sum)
}

#' @title Print PS estimation results
#'
#' @description Print summary information of PS estimation results
#'
#' @param x A list of class \code{PSRWE_DTA_MAT} that is generated using
#'   the \code{\link{psrwe_match}} function.
#' @param ... Additional parameters
#'
#' @seealso  \code{\link{summary.PSRWE_DTA_MAT}}
#'
#'
#' @method print PSRWE_DTA_MAT
#'
#' @return A list from \code{summary(x)} with additional information
#'
#' @export
#'
print.PSRWE_DTA_MAT <- function(x, ...) {
    rst_sum <- summary(x, ...)

    ## overall summary
    cat_ps_dta(x, rst_sum)
    cat("\n")

    ## matched
    ss <- NULL
    if (!is.null(x$strata_covs)) {
        ss <- paste("The matching is stratified by the following covariate(s): ",
                    paste(x$strata_covs, collapse = ", "),
                    ". ", sep = "")
    }

    ss <- paste(ss, "A total of ",
                sum(rst_sum$N_Match[c("Zero", "Other")]),
                " current study subjects are matched by less than ",
                rst_sum$ratio,
                " RWD subjects. Please note unequal matching may cause",
                " unbalance in covariate distributions of the current",
                " and matched RWD subjects.", sep = "")

    cat_paste(ss)
    cat("\n")

    ## summary table
    cat_paste("The following table summarizes the number of subjects in each stratum:")
    cat("\n")

    print(rst_sum$Summary)
}


#' Plot PS distributions
#'
#' S3 method for visualizing PS adjustment based on matching.
#'
#' @param x A list of class \code{PSRWE_DTA_MAT} that is generated using
#'   the \code{\link{psrwe_match}} function.
#'
#' @param ... Parameters for \code{plot.PSRWE_DAT}
#'
#' @seealso  \code{\link{plot.PSRWE_DTA}}
#'
#' @method plot PSRWE_DTA_MAT
#'
#' @return A plot of class in ggplot2
#'
#' @export
#'
plot.PSRWE_DTA_MAT <- function(x, ...) {
    plot.PSRWE_DTA(x, ...)
}
