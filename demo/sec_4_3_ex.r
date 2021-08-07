### Example of Section 4.3.
suppressMessages(library(psrwe, quietly = TRUE))
data(ex_dta)

### First parts of Data.
head(ex_dta)

### Obtain PSs.
dta_ps_single <- rwe_ps_est(ex_dta,
                     v_covs = paste("V", 1:7, sep = ""),
                     v_grp = "Group", cur_grp_level = "current",
                     ps_method = "logistic", nstrata = 5)

### Sample size.
### See "sec_4_1_ex" for details.

### Balance assessment of PS stratification.
### See "sec_4_1_ex" for details.

### Obtain discounting parameters.
### See "sec_4_1_ex" for details.
ps_bor_single <- rwe_ps_borrow(dta_ps_single, total_borrow = 30,
                               method = "distance", metric = "ovl")

### PSKM, single arm study, time-to-event outcome.
rst_km <- rwe_ps_survkm(ps_bor_single,
                        v_time    = "Y_Surv",
                        v_event   = "Status",
                        pred_tp  = 365)
rst_km

### Plot PSKM
rst_km_allt <- rwe_ps_survkmplot(ps_bor_single,
                                 v_time    = "Y_Surv",
                                 v_event   = "Status")
rst_km_ci <- rwe_ps_survkmci(rst_km_allt)
plot_survkm(rst_km_ci)
