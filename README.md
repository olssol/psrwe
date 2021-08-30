
<!-- README.md is generated from README.Rmd. Please edit that file -->

# psrwe

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/psrwe)](https://CRAN.R-project.org/package=psrwe)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![Download](https://cranlogs.r-pkg.org/badges/psrwe)](https://cran.r-project.org/package=psrwe)
[![Build
Status](https://app.travis-ci.com/snoweye/psrwe.svg?branch=ci)](https://app.travis-ci.com/snoweye/psrwe)
[![Appveyor
Build status](https://ci.appveyor.com/api/projects/status/lnta05dn3ex9v641?svg=true)](https://ci.appveyor.com/project/snoweye/psrwe)
<!-- badges: end -->

High-quality real-world data can be transformed into scientific
real-world evidence (RWE) for regulatory and healthcare decision-making
using proven analytical methods and techniques. For example, propensity
score (PS) methodology can be applied to pre-select a subset of
real-world data containing patients that are similar to those in the
current clinical study in terms of covariates, and to stratify the
selected patients together with those in the current study into more
homogeneous strata. Then, methods such as the power prior approach or
composite likelihood approach can be applied in each stratum to draw
inference for the parameters of interest. This package provides
functions that implement the PS-integrated RWE analysis methods proposed
in [Wang et al. (2019)](https://doi.org/10.1080/10543406.2019.1657133),
[Wang et al. (2020)](https://doi.org/10.1080/10543406.2019.1684309), and
[Chen et al. (2020)](https://doi.org/10.1080/10543406.2020.1730877).

## Installation

You can install the released version of `psrwe` from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("psrwe")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("olssol/psrwe")
```

## References

1.  Wang C, Li H, Chen WC, Lu N, Tiwari R, Xu Y, Yue LQ. Propensity
    score-integrated power prior approach for incorporating real-world
    evidence in single-arm clinical studies. *Journal of
    Biopharmaceutical Statistics*, 2019; **29**, 731–748.
    <https://doi.org/10.1080/10543406.2019.1657133>.

2.  Chen WC, Wang C, Li H, Lu N, Tiwari R, Xu Y, Yue LQ. (2020),
    Propensity score-integrated composite likelihood approach for
    augmenting the control arm of a randomized controlled trial by
    incorporating real-world data. *Journal of Biopharmaceutical
    Statistics*, 2020; **30**, 508–520.
    <https://doi.org/10.1080/10543406.2020.1730877>.

3.  Wang C, Lu N, Chen WC, Li H, Tiwari R, Xu Y, Yue LQ. (2020),
    Propensity score-integrated composite likelihood approach for
    incorporating real-world evidence in single-arm clinical studies.
    *Journal of Biopharmaceutical Statistics*, 2020; **30**, 495–507.
    <https://doi.org/10.1080/10543406.2019.1684309>.
