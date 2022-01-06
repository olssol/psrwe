#' PS-Integrated Methods for Incorporating RWE in Clinical Studies
#'
#' @docType   package
#' @name      psrwe-package
#' @aliases   psrwe
#' @useDynLib psrwe, .registration = TRUE
#'
#' @import Rcpp
#' @import methods
#' @import ggplot2
#' @import rstantools
#' @importFrom stats approxfun as.formula binomial cov density ecdf glm
#'     integrate optim predict quantile sd var ks.test qnorm pnorm
#' @importFrom rstan sampling extract stanc rstan_options traceplot stan_rhat
#' @importFrom randomForest randomForest
#' @importFrom grDevices colors
#' @importFrom graphics axis box legend lines par plot points text
#' @importFrom parallel detectCores
#' @importFrom cowplot plot_grid
#' @importFrom dplyr %>% group_by_ group_by summarize mutate count mutate_if
#'     rename filter select arrange ungroup n distinct left_join if_else
#' @importFrom survival Surv survfit
#'
#' @description
#'
#' This package provide R functions for conducting clinical studies with
#' real-world evidence (RWE) incorporated in the study design and analysis.
#'
#' @section PS-integrated power prior:
#'
#' We extend the Bayesian power prior approach for a single-arm study (the
#' current study) to leverage external real-world data (RWD). We use propensity
#' score methodology to pre-select a subset of real-world data containing
#' patients that are similar to those in the current study in terms of
#' covariates, and to stratify the selected patients together with those in the
#' current study into more homogeneous strata. The power prior approach is then
#' applied in each stratum to obtain stratum-specific posterior distributions,
#' which are combined to complete the Bayesian inference for the parameters of
#' interest.
#'
#' @section PS-integrated composite likelihood:
#'
#' A propensity score-integrated composite likelihood (PSCL) approach is
#' developed for cases in which the control arm of a two-arm randomized
#' controlled trial (RCT) (treated vs. control) is augmented with patients from
#' real-world data (RWD) containing both clinical outcomes and covariates at the
#' patient-level. The PSCL approach first estimates the propensity score for
#' every patient as the probability of the patient being in the RCT rather than
#' the RWD, and then stratifies all patients into strata based on the estimated
#' propensity scores. Within each propensity score stratum, a composite
#' likelihood function is specified and utilized to down-weight the information
#' contributed by the RWD source. Estimates of the stratum-specific parameters
#' are obtained by maximizing the composite likelihood function. These
#' stratum-specific estimates are then combined to obtain an overall
#' population-level estimate of the parameter of interest.
#'
#' @references
#'
#' Chen WC, Wang C, Li H, Lu N, Tiwari R, Xu Y, Yue LQ. Propensity
#' score-integrated composite likelihood approach for augmenting the control arm
#' of a randomized controlled trial by incorporating real-world data.
#' *Journal of Biopharmaceutical Statistics*. 2020; 30(3):508-520.
#'
#' Wang C, Lu N, Chen WC, Li H, Tiwari R, Xu Y, Yue LQ. Propensity
#' score-integrated composite likelihood approach for incorporating real-world
#' evidence in single-arm clinical studies.
#' *Journal of Biopharmaceutical Statistics*. 2020; 30(3):495-507.
#'
#' Wang C, Li H, Chen WC, Lu N, Tiwari R, Xu Y, Yue LQ. Propensity
#' score-integrated power prior approach for incorporating real-world evidence
#' in single-arm clinical studies. *Journal of Biopharmaceutical Statistics*.
#' 2019; 29(5):731-748.
NULL


#' @title Example dataset
#'
#' Example dataset of a single arm study.
#'
#' @usage data(ex_dta)
#' @name ex_dta
#' @keywords datasets
#'
#' @format A data frame with the following variables:
#' \itemize{
#'   \item{Group}{current, rwd}
#'   \item{Y_Bin}{Binary outcome}
#'   \item{Y_Con}{Continuous outcome}
#'   \item{Y_Surv}{Survival outcome in days}
#'   \item{Status}{Event status (0=alive, 1=dead)}
#'   \item{V1-V7}{Covariates}
#' }
"ex_dta"


#' @title Example dataset
#'
#' Example dataset of a randomized study.
#'
#' @usage data(ex_dta_rct)
#' @name ex_dta_rct
#' @keywords datasets
#'
#' @format A data frame with the following variables:
#' \itemize{
#'   \item{Group}{current, rwd}
#'   \item{Arm}{control, treatment}
#'   \item{Y_Con}{Continuous outcome}
#'   \item{V1-V7}{Covariates}
#' }
"ex_dta_rct"
