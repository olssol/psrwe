### Example of RCT, time-to-event outcome, and log-rank test
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
rst_lrk <- psrwe_survlrk(ps_bor_rct,
                         pred_tp = 365,
                         v_time = "Y_Surv",
                         v_event = "Status")

### Outcome analysis.
oa_lrk <- psrwe_outana(rst_lrk)
oa_lrk
print(oa_lrk, show_details = TRUE)
summary(oa_lrk, pred_tps = c(180, 365))

### Use simple Bootstrap stderr. This may take a while longer.
set.seed(12341)
rst_lrk_sbs <- psrwe_survlrk(ps_bor_rct,
                             pred_tp = 365,
                             v_time = "Y_Surv",
                             v_event = "Status",
                             stderr_method = "sbs")
oa_lrk_sbs <- psrwe_outana(rst_lrk_sbs)
summary(oa_lrk_sbs, pred_tps = c(180, 365))

### Use complex Bootstrap stderr. This may take a while longer.
set.seed(12342)
rst_lrk_cbs <- psrwe_survlrk(ps_bor_rct,
                             pred_tp = 365,
                             v_time = "Y_Surv",
                             v_event = "Status",
                             stderr_method = "cbs")
oa_lrk_cbs <- psrwe_outana(rst_lrk_cbs)
summary(oa_lrk_cbs, pred_tps = c(180, 365))

## Reset to user's options.
options(org_digits)

