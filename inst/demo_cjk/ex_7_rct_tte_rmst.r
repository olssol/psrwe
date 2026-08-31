### Example of RCT, time-to-event outcome, and RMST
suppressMessages(library(psrwe, quietly = TRUE))
org_digits <- options(digits = 3)
data(ex_dta_rct)

### Obtain PSs.
dta_ps_rct <- psrwe_est(ex_dta_rct,
                        v_covs = paste("V", 1:7, sep = ""),
                        v_grp = "Group", cur_grp_level = "current",
                        v_arm = "Arm", ctl_arm_level = "control",
                        ps_method = "logistic", nstrata = 5,
                        stra_ctl_only = FALSE)

### Obtain discounting parameters.
ps_bor_rct <- psrwe_borrow(dta_ps_rct, total_borrow = 30)

### PSLRK, two-arm RCT, time-to-event outcome.
rst_rmst <- psrwe_survrmst(ps_bor_rct,
                           pred_tp = 365,
                           v_time = "Y_Surv",
                           v_event = "Status")

### Outcome analysis.
oa_rmst <- psrwe_outana(rst_rmst)
oa_rmst
print(oa_rmst, show_details = TRUE)
summary(oa_rmst, pred_tps = c(180, 365))

### Use complex Jackknife stderr. This may take a while longer.
rst_rmst_cjk <- psrwe_survrmst(ps_bor_rct,
                               pred_tp = 365,
                               v_time = "Y_Surv",
                               v_event = "Status",
                               stderr_method = "cjk")
oa_rmst_cjk <- psrwe_outana(rst_rmst_cjk)
summary(oa_rmst_cjk, pred_tps = c(180, 365))

## Reset to user's options.
options(org_digits)

