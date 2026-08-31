#' @title Distance between two distributions
#'
#' @description Calculate difference measures using different metrics.
#'
#' @param cov0 Vector (or matrix for \code{metric = "mhb"}) of samples from the
#'   first distribution.
#' @param cov1 Vector (or matrix for \code{metric = "mhb"}) of samples from the
#'   second distribution.
#' @param metric Metric to use for calculating the distance with options:
#' \describe{
#'   \item{\code{ovl}}{Overlapping area} (default)
#'   \item{\code{ksd}}{Kullback-Leibler distance}
#'   \item{\code{astd}}{Standardized absolute mean difference}
#'   \item{\code{std}}{Standardized mean difference}
#'   \item{\code{abd}}{Absolute difference in means}
#'   \item{\code{ley}}{Levy distance}
#'   \item{\code{mhb}}{Mahalanobis distance}
#'   \item{\code{omkss}}{One minus Kolmogorov-Smirnov statistic}
#' }
#'
#' @return A real value of the distance.
#'
#' @examples
#'
#' x <- rnorm(100,  mean = 0, sd = 1)
#' y <- rnorm(1000, mean = 1, sd = 2)
#' get_distance(x, y, "ovl")
#' get_distance(x, y, "abd")
#'
#' @export
#'
get_distance <- function(cov0, cov1,
                         metric = c("ovl", "ksd", "astd", "std", "abd",
                                    "ley", "mhb", "omkss")) {
    metric <- match.arg(metric)
    switch(metric,
           abd   = abs(mean(cov0) - mean(cov1)),
           std   = metric_std(cov0, cov1),
           astd  = metric_std(cov0, cov1, TRUE),
           ovl   = metric_ovl(cov0, cov1),
           ksd   = metric_ksd(cov0, cov1),
           ley   = metric_ley(cov0, cov1),
           mhb   = metric_mhb(cov0, cov1),
           omkss = metric_omkss(cov0, cov1)
           )
}


#' @title Overlapping coefficient
#'
#' @inheritParams get_distance
#'
#' @noRd
metric_ovl <- function(cov0, cov1) {
  cov <- c(cov0, cov1)
  if (length(unique(cov)) <= 10) {
      all_x <- c(rep(0, length(cov0)),
                 rep(1, length(cov1)))
      pt    <- apply(prop.table(table(cov, all_x), 2), 1, min)

    return(sum(pt))
  }

  mn <- pmax(0, min(cov) - 1e-3)
  mx <- pmin(1, max(cov) + 1e-3)

  f1 <- approxfun(density(cov1, from = mn, to = mx,
                          bw = "nrd"))
  f0 <- approxfun(density(cov0, from = mn, to = mx,
                          bw = "nrd"))

  s <- try(integrate(function(x) pmin(f1(x), f0(x)),
                     lower = mn, upper = mx,
                     subdivisions = 500)$value)

  ifelse(inherits(s, "try-error"), NA, s)
}

#' @title K-S distance
#'
#' @inheritParams get_distance
#'
#' @noRd
metric_ksd <- function(cov0, cov1) {
    cov    <- c(cov0, cov1)
    cdf_1  <- ecdf(cov1)
    cdf_0  <- ecdf(cov0)
    1 / max(abs(cdf_1(cov) - cdf_0(cov)))
}

#' @title Levy distance
#'
#' @inheritParams get_distance
#'
#' @noRd
metric_ley <- function(cov0, cov1) {
  cov   <- c(cov0, cov1)
  cdf_1 <- ecdf(cov1)
  cdf_0 <- ecdf(cov0)
  e     <- max(abs(cdf_1(cov) - cdf_0(cov)))

  if (length(unique(cov)) <= 10) {
    return(e)
  }

  x     <- seq(min(cov), max(cov), length.out = 1000)
  check <- all(cdf_0(x - e) - e <= cdf_1(x) &
                 cdf_1(x) <= cdf_0(x + e) + e)

  while (check) {
    e <- e - 1e-12
    check <- all(cdf_0(x - e) - e <= cdf_1(x) &
                   cdf_1(x) <= cdf_0(x + e) + e)
  }

  1 / e
}

#' @title Mahalanobis balance
#'
#' @inheritParams get_distance
#'
#' @details Both \code{covs} should be a reduced dataset that contains only those
#'   covariates that will be used for calculating Mahalanobis balance, for
#'   example, \code{covs = dat[, 1:6]}. \code{trt} should be the exposure
#'   variable, for example, \code{trt = dat$X}.
#'
#' @noRd
metric_mhb <- function(cov0, cov1) {
  cov01 <- rbind(cov0, cov1)
  sinv  <- solve(cov(cov01))
  x0    <- colMeans(cov0)
  x1    <- colMeans(cov1)

  rst <- sum((t(x1 - x0) %*% sinv) * (x1 - x0))
  1 / rst
}

#' @title One minus Kolmogorov-Smirnov statistic
#'
#' @inheritParams get_distance
#'
#' @noRd
metric_omkss <- function(cov0, cov1) {
    1 - ks.test(cov0[!is.na(cov0)],
                cov1[!is.na(cov0)])$statistic
}

#' @title Standardized difference in means
#'
#' @inheritParams get_distance
#'
#' @param is_abs Whether to take the absolute values
#'
#' @noRd
metric_std <- function(cov0, cov1, is_abs = FALSE) {
    s   <- sqrt((var(cov1) + var(cov0)) / 2)
    rst <- (mean(cov1) - mean(cov0)) / s

    if (is_abs) {
        rst <- abs(rst)
    }

    rst
}
