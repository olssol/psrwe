## ----eval=T, echo=FALSE-------------------------------------------------------
suppressMessages(require(psrwe, quietly = TRUE))
org_digits <- options(digits = 3)
set.seed(1000)

## ----eval=T, echo=TRUE--------------------------------------------------------
data(ex_dta)
dta_ps <- psrwe_est(ex_dta,
                    v_covs = paste("V", 1:7, sep = ""),
                    v_grp = "Group",
                    cur_grp_level = "current",
                    ps_method = "logistic")
dta_ps

## ----eval=T, echo=TRUE--------------------------------------------------------
dta_ps_match <- psrwe_match(dta_ps,
                            ratio = 2,
                            strata_covs = "V1",
                            seed = 123)
dta_ps_match

## ----eval=T, echo=TRUE--------------------------------------------------------
ps_bor_match <- psrwe_borrow(dta_ps_match,
                             total_borrow = 30)
ps_bor_match

## ----eval=T, echo=TRUE--------------------------------------------------------
rst_cl <- psrwe_compl(ps_bor_match,
                      outcome_type = "binary",
                      v_outcome    = "Y_Bin")
rst_cl

## ----eval=T, echo=TRUE--------------------------------------------------------
oa_cl <- psrwe_outana(rst_cl, method_ci = "wilson", mu = 0.40)
oa_cl

## ----eval=T, echo=FALSE-------------------------------------------------------
## Reset to user's options.
options(org_digits)

