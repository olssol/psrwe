### Example of Section 4.4.
suppressMessages(library(psrwe, quietly = TRUE))
data(ex_dta)

### First parts of Data.
head(ex_dta)

### Obtain PSs.
dta_ps_single <- ps_rwe_est(ex_dta,
                     v_covs = paste("V", 1:7, sep = ""),
                     v_grp = "Group", cur_grp_level = "current",
                     ps_method = "logistic", nstrata = 5)

### Balance assessment of PS stratification.
### See "sec_4_1_ex" for details.

### Obtain discounting parameters.
### See "sec_4_1_ex" for details.
ps_bor_single <- ps_rwe_borrow(dta_ps_single, total_borrow = 30)

### PSKM, single arm study, time-to-event outcome.
rst_km <- ps_rwe_survkm(ps_bor_single,
                        v_time    = "Y_Surv",
                        v_event   = "Status",
                        pred_tp  = 365)
rst_km

### Plot PSKM
plot(rst_km)
plot(rst_km, add_ci = FALSE, add_stratum = TRUE)
plot(rst_km, conf_type = "plain")

### 95% two-sided CI
rst_km_log <- ps_rwe_ci(rst_km)
rst_km_log
rst_km_pln <- ps_rwe_ci(rst_km, conf_type = "plain")
rst_km_pln

### Inference.
rst_km <- ps_rwe_infer(rst_km, mu = 0.70, alternative = "greater")
rst_km

