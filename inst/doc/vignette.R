## ----eval=T, echo=FALSE-------------------------------------------------------
require(psrwe)
set.seed(1000)

## ----eval=T, echo=TRUE--------------------------------------------------------
data(ex_dta)
dta_ps <- psrwe_est(ex_dta,
                     v_covs = paste("V", 1:7, sep = ""),
                     v_grp = "Group",
                     cur_grp_level = "current",
                     nstrata = 5,
                     ps_method = "logistic")
dta_ps

## ----echo=TRUE, fig.width=6, fig.height=5-------------------------------------
plot(dta_ps, plot_type = "balance")

## ----echo=TRUE, fig.width=6, fig.height=5-------------------------------------
plot(dta_ps, plot_type = "ps")

## ----eval=T, echo=TRUE--------------------------------------------------------
ps_bor <- psrwe_borrow(dta_ps, total_borrow = 40,
                        method = "distance")
rst_pp <- psrwe_powerp(ps_bor, v_outcome = "Y_Bin",
                        outcome_type = "binary",
                        seed = 1234)

## ----eval=T, echo=TRUE--------------------------------------------------------
summary(rst_pp)

## ----eval=T, echo=TRUE--------------------------------------------------------
rst_cl <- psrwe_compl(ps_bor, v_outcome = "Y_Bin",
                       outcome_type = "binary")
summary(rst_cl)

## ----eval=T, echo=TRUE--------------------------------------------------------
data(ex_dta_rct)
dta_ps_rct <- psrwe_est(ex_dta_rct, v_covs = paste("V", 1:7, sep = ""),
                         v_grp = "Group", cur_grp_level = "current",
                         v_arm = "Arm", ctl_arm_level = "control")
dta_ps_rct

ps_bor_rct <- psrwe_borrow(dta_ps_rct, total_borrow = 30,
                            method = "distance")
ps_bor_rct

rst_cl_rct <- psrwe_compl(ps_bor_rct, v_outcome = "Y_Con",
                           outcome_type = "continuous")

rst_cl_rct$Effect

