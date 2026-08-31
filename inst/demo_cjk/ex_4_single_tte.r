### Example of single arm and time-to-event outcome
suppressMessages(library(psrwe, quietly = TRUE))
org_digits <- options(digits = 3)
data(ex_dta)

### Obtain PSs.
dta_ps_single <- psrwe_est(ex_dta,
                     v_covs = paste("V", 1:7, sep = ""),
                     v_grp = "Group", cur_grp_level = "current",
                     ps_method = "logistic", nstrata = 5)

### Obtain discounting parameters.
ps_bor_single <- psrwe_borrow(dta_ps_single, total_borrow = 30)

### PSKM, single arm study, time-to-event outcome.
rst_km <- psrwe_survkm(ps_bor_single,
                       pred_tp  = 365,
                       v_time    = "Y_Surv",
                       v_event   = "Status")

### Outcome analysis.
oa_km <- psrwe_outana(rst_km, mu = 0.70, alternative = "greater")
oa_km
print(oa_km, show_details = TRUE)
summary(oa_km, pred_tps = c(180, 365))

### Use complex Jackknife stderr. This may take a while longer.
rst_km_cjk <- psrwe_survkm(ps_bor_single,
                           pred_tp  = 365,
                           v_time    = "Y_Surv",
                           v_event   = "Status",
                           stderr_method = "cjk")
oa_km_cjk <- psrwe_outana(rst_km_cjk, mu = 0.70, alternative = "greater")
summary(oa_km_cjk, pred_tps = c(180, 365))

## Reset to user's options.
options(org_digits)

