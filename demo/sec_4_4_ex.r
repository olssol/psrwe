### Example of Section 4.4.
suppressMessages(library(psrwe, quietly = TRUE))
options(digits = 3)
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
                       v_time    = "Y_Surv",
                       v_event   = "Status",
                       pred_tp  = 365)
rst_km

### Plot PSKM
plot(rst_km)
plot(rst_km, add_ci = FALSE, add_stratum = TRUE)
plot(rst_km, conf_type = "plain")

### Outcome analysis.
oa_km <- psrwe_outana(rst_km, mu = 0.70, alternative = "greater")
oa_km

