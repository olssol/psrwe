### Example of Section 4.5.
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

### PSKM, two-arm RCT, time-to-event outcome.
rst_km_rct <- psrwe_survkm(ps_bor_rct,
                           pred_tp = 365,
                           v_time = "Y_Surv",
                           v_event = "Status")
rst_km_rct

### Plot PSKM.
plot(rst_km_rct, xlim = c(0, 730))

### Outcome analysis.
oa_km_rct <- psrwe_outana(rst_km_rct, alternative = "greater")
oa_km_rct
print(oa_km_rct, show_rct = TRUE)
summary(oa_km_rct, pred_tps = c(180, 365))

### Use Jackknife stderr. This may take a while.
rst_km_rct_jk <- psrwe_survkm(ps_bor_rct,
                              pred_tp = 365,
                              v_time = "Y_Surv",
                              v_event = "Status",
                              stderr_method = "jk")
oa_km_rct_jk <- psrwe_outana(rst_km_rct_jk, alternative = "greater")
summary(oa_km_rct_jk, pred_tps = c(180, 365))

### Use simple Jackknife stderr. This may take a while longer.
rst_km_rct_jko <- psrwe_survkm(ps_bor_rct,
                               pred_tp = 365,
                               v_time = "Y_Surv",
                               v_event = "Status",
                               stderr_method = "sjk")
oa_km_rct_jko <- psrwe_outana(rst_km_rct_jko, alternative = "greater")
summary(oa_km_rct_jko, pred_tps = c(180, 365))

### Reset to user's options.
options(org_digits)

