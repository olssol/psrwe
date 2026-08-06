## ----eval=T, echo=FALSE-------------------------------------------------------
suppressMessages(require(psrwe, quietly = TRUE))
org_digits <- options(digits = 3)
set.seed(1000)

## ----eval=T, echo=TRUE--------------------------------------------------------
data(ex_dta_rct)
dta_ps_rct <- psrwe_est(ex_dta_rct,
                        v_covs = paste("V", 1:7, sep = ""),
                        v_grp = "Group", cur_grp_level = "current",
                        v_arm = "Arm", ctl_arm_level = "control",
                        ps_method = "logistic", nstrata = 5,
                        stra_ctl_only = FALSE)
ps_bor_rct <- psrwe_borrow(dta_ps_rct, total_borrow = 30)

## ----eval=T, echo=TRUE--------------------------------------------------------
rst_km_rct <- psrwe_survkm(ps_bor_rct,
                           pred_tp = 365,
                           v_time = "Y_Surv",
                           v_event = "Status")
rst_km_rct

## ----echo=TRUE, fig.width=6, fig.height=5-------------------------------------
plot(rst_km_rct, xlim = c(0, 730))

## ----eval=T, echo=TRUE--------------------------------------------------------
oa_km_rct <- psrwe_outana(rst_km_rct, alternative = "greater")
oa_km_rct

## ----eval=T, echo=TRUE--------------------------------------------------------
print(oa_km_rct, show_rct = TRUE)

## ----eval=T, echo=TRUE--------------------------------------------------------
summary(oa_km_rct, pred_tps = c(180, 365))

## ----eval=T, echo=TRUE--------------------------------------------------------
rst_lrk <- psrwe_survlrk(ps_bor_rct,
                         pred_tp = 365,
                         v_time = "Y_Surv",
                         v_event = "Status")
rst_lrk

## ----eval=T, echo=TRUE--------------------------------------------------------
oa_lrk <- psrwe_outana(rst_lrk)
oa_lrk

## ----eval=T, echo=TRUE--------------------------------------------------------
print(oa_lrk, show_details = TRUE)

## ----eval=T, echo=TRUE--------------------------------------------------------
summary(oa_lrk, pred_tps = c(180, 365))

## ----eval=T, echo=TRUE--------------------------------------------------------
rst_rmst <- psrwe_survrmst(ps_bor_rct,
                           pred_tp = 365,
                           v_time = "Y_Surv",
                           v_event = "Status")
rst_rmst

## ----eval=T, echo=TRUE--------------------------------------------------------
oa_rmst <- psrwe_outana(rst_rmst)
oa_rmst

## ----eval=T, echo=TRUE--------------------------------------------------------
print(oa_rmst, show_details = TRUE)

## ----eval=T, echo=TRUE--------------------------------------------------------
summary(oa_rmst, pred_tps = c(180, 365))

## ----eval=T, echo=FALSE-------------------------------------------------------
## Reset to user's options.
options(org_digits)

