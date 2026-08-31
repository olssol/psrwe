### Example of Section 4.4.
suppressMessages(library(psrwe, quietly = TRUE))
org_options <- options(digits = 3)
data(ex_dta)

### First parts of Data.
head(ex_dta)

### Obtain PSs.
dta_ps_single <- psrwe_est(ex_dta,
                     v_covs = paste("V", 1:7, sep = ""),
                     v_grp = "Group", cur_grp_level = "current",
                     ps_method = "logistic", nstrata = 5)

### Balance assessment of PS stratification.
### See "sec_4_1_ex" for details.

### Obtain discounting parameters.
### See "sec_4_1_ex" for details.
ps_bor_single <- psrwe_borrow(dta_ps_single, total_borrow = 30)

### PSKM, single arm study, time-to-event outcome.
rst_km <- psrwe_survkm(ps_bor_single,
                       pred_tp  = 365,
                       v_time    = "Y_Surv",
                       v_event   = "Status")
rst_km

### Plot PSKM.
plot(rst_km)
plot(rst_km, add_ci = FALSE, add_stratum = TRUE, ylim = c(0.35, 1))
plot(rst_km, conf_type = "plain", ylim = c(0.35, 1))

### Outcome analysis.
oa_km <- psrwe_outana(rst_km, mu = 0.70, alternative = "greater")
oa_km
print(oa_km, show_details = TRUE)
summary(oa_km, pred_tps = c(180, 365))

### Use Jackknife stderr. This may take a while.
rst_km_jk <- psrwe_survkm(ps_bor_single,
                          pred_tp  = 365,
                          v_time    = "Y_Surv",
                          v_event   = "Status",
                          stderr_method = "jk")
oa_km_jk <- psrwe_outana(rst_km_jk, mu = 0.70, alternative = "greater")
summary(oa_km_jk, pred_tps = c(180, 365))

### Use simple Jackknife stderr. This may take a while longer.
rst_km_jko <- psrwe_survkm(ps_bor_single,
                           pred_tp  = 365,
                           v_time    = "Y_Surv",
                           v_event   = "Status",
                           stderr_method = "sjk")
oa_km_jko <- psrwe_outana(rst_km_jko, mu = 0.70, alternative = "greater")
summary(oa_km_jko, pred_tps = c(180, 365))

### Reset to user's options.
options(org_digits)

