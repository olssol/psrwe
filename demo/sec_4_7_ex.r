### Example of Section 4.7.
suppressMessages(library(psrwe, quietly = TRUE))
org_options <- options(digits = 3)
data(ex_dta_rct)

### First parts of Data.
head(ex_dta_rct)

### Obtain PSs.
dta_ps_rct <- psrwe_est(ex_dta_rct,
                        v_covs = paste("V", 1:7, sep = ""),
                        v_grp = "Group", cur_grp_level = "current",
                        v_arm = "Arm", ctl_arm_level = "control",
                        ps_method = "logistic", nstrata = 5,
                        stra_ctl_only = FALSE)

### Balance assessment of PS stratification.
### See "sec_4_2_ex" for details.

### Obtain discounting parameters.
### See "sec_4_2_ex" for details.
ps_bor_rct <- psrwe_borrow(dta_ps_rct, total_borrow = 30)

### PSLRK, two-arm RCT, time-to-event outcome.
rst_rmst <- psrwe_survrmst(ps_bor_rct,
                           pred_tp = 365,
                           v_time = "Y_Surv",
                           v_event = "Status")
rst_rmst

### Outcome analysis.
oa_rmst <- psrwe_outana(rst_rmst)
oa_rmst
print(oa_rmst, show_details = TRUE)
summary(oa_rmst, pred_tps = c(180, 365))

### Use Jackknife stderr. This may take a while.
rst_rmst_jk <- psrwe_survrmst(ps_bor_rct,
                              pred_tp = 365,
                              v_time = "Y_Surv",
                              v_event = "Status",
                              stderr_method = "jk")
oa_rmst_jk <- psrwe_outana(rst_rmst_jk)
summary(oa_rmst_jk, pred_tps = c(180, 365))

### Use simple Jackknife stderr. This may take a while longer.
rst_rmst_jko <- psrwe_survrmst(ps_bor_rct,
                               pred_tp = 365,
                               v_time = "Y_Surv",
                               v_event = "Status",
                               stderr_method = "sjk")
oa_rmst_jko <- psrwe_outana(rst_rmst_jko)
summary(oa_rmst_jko, pred_tps = c(180, 365))

### Reset to user's options.
options(org_digits)

