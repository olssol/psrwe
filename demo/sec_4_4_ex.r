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

### Balance assessment of PS stratification.
### See "sec_4_1_ex" for details.

### Obtain discounting parameters.
### See "sec_4_1_ex" for details.
ps_bor_single <- rwe_ps_borrow(dta_ps_single, total_borrow = 30)

### PSKM, single arm study, time-to-event outcome.
rst_km <- rwe_ps_survkm(ps_bor_single,
                        v_time    = "Y_Surv",
                        v_event   = "Status",
                        pred_tp  = 365)
rst_km

### Plot PSKM
plot(rst_km)
