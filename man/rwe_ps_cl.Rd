% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/psrwe_complikel.R
\name{rwe_ps_cl}
\alias{rwe_ps_cl}
\title{PS-Integrated Composite Likelihood Estimation}
\usage{
rwe_ps_cl(data_withps, v_borrow = 0, v_outcome = "Y", ...)
}
\arguments{
\item{data_withps}{A list of class \code{RWE_DWITHPS} that is generated using
the \code{\link{rwe_ps}} function.}

\item{v_borrow}{Vector of number of subjects to be borrowed from each stratum}

\item{v_outcome}{Column name corresponding to the outcome}

\item{...}{Parameters for \code{rwe_cl}}
}
\value{
A dataframe with class name \code{RWE_CLRST}. It contains the
    composite estimation of the mean for each stratum as well as the
    jackknife estimation for each subject. The results should be further
    summarized by its S3 method \code{summary}.
}
\description{
Estimate the mean of the outcome based on PS-integrated composite likelihood
approach. Variance is estimated by Jack-Knife method. Applies to the case
when there is only one external data source.
}
\examples{
\donttest{
dta_ps <- rwe_ps(ex_dta,
                 v_covs = paste("V", 1:7, sep = ""),
                 v_grp = "Group",
                 cur_grp_level = "current",
                 nstrata = 5)
ps_dist <- rwe_ps_dist(dta_ps)
ps_borrow <- rwe_ps_borrow(total_borrow = 40, ps_dist)
rst_cl    <- rwe_ps_cl(dta_ps, v_borrow = ps_borrow)
summary(rst_cl)}

}
