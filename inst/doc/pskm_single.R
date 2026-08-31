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
                    nstrata = 5,
                    ps_method = "logistic")
ps_bor <- psrwe_borrow(dta_ps,
                       total_borrow = 30,
                       method = "distance")

## ----eval=T, echo=TRUE--------------------------------------------------------
rst_km <- psrwe_survkm(ps_bor,
                       pred_tp = 365,
                       v_time  = "Y_Surv",
                       v_event = "Status")
rst_km

## ----echo=TRUE, fig.width=6, fig.height=5-------------------------------------
plot(rst_km)

## ----echo=TRUE, fig.width=6, fig.height=5-------------------------------------
plot(rst_km, add_ci = FALSE, add_stratum = TRUE)

## ----echo=TRUE, fig.width=6, fig.height=5-------------------------------------
plot(rst_km, conf_type = "plain")

## ----eval=T, echo=TRUE--------------------------------------------------------
oa_km <- psrwe_outana(rst_km, mu = 0.70, alternative = "greater")
oa_km

## ----eval=T, echo=TRUE--------------------------------------------------------
print(oa_km, show_details = TRUE)

## ----eval=T, echo=TRUE--------------------------------------------------------
summary(oa_km, pred_tps = c(180, 365))

## ----eval=T, echo=FALSE-------------------------------------------------------
## Reset to user's options.
options(org_digits)

