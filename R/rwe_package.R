#' R tools for synthesizing real-world evidence
#'
#' @docType   package
#' @name      psrwe-package
#' @aliases   psrwe
#' @useDynLib psrwe, .registration = TRUE
#'
#' @import stats
#' @import Rcpp
#' @import methods
#' @import ggplot2
#'
#' @importFrom rstan         sampling extract stanc rstan_options traceplot stan_rhat
#' @importFrom randomForest  randomForest
#'
#' @importFrom grDevices colors
#' @importFrom graphics  axis box legend lines par plot points text
#'
#' @importFrom parallel  detectCores
#' @importFrom mvtnorm   rmvnorm
#' @importFrom cowplot   plot_grid
#' @importFrom dplyr     %>% group_by_ group_by summarize mutate count_ mutate_if rename_ filter
#'
#'
#' @description
#'
#' This package contains the functions for synthesizing real-world evidence in
#' single-arm medical device studies.
NULL



#' Parameters for simulating data
#'
#' @name simupara
#'
#' @param nPat number of patients
#' @param muCov mean vector of covariates
#' @param sdCov standard deviation vector of covariates
#' @param corCov correlation of covariates
#' @param regCoeff regression coefficients
#' @param mix.phi weight in mixture model
#' @param cov.breaks breaks to cut covaraites into categories. If numeric, the
#'     same cuts will be applied to all covariates. If list, each covariates
#'     will be categorized based on its own breaks. If NULL, keep continuous.
#' @param type distributions of the random error
#' @param ysig standard error of the random error
#' @param skew.n parameter of negative bionomial distribution
#' @param skew.p parameter of negative bionomial distribution
#' @param b0 intercept in regession model
#' @param bin.mu mean of the binary outcomes used to compute b0
#' @param formula.z formula of the treatement assignment model. No intercept
#'     term.
#' @param formula.y formula of the outcome model. No intercept term.
#' @param trial.data existing clinical trial data
#' @param group column referring to arm in the existing dataset
#' @param outcome column referring to outcome in the existing dataset
#' @param trt.effec true treatment effect to be set
#' @param with.replacement sample with or without replacement from the existing
#'     dataset
#' @param keep.group sample ignore arm assignment in the existing dataset
#'
#'
NULL
