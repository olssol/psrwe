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
table(dta_ps_single$data$"_grp_")     # by data source
table(dta_ps_single$data$"_strata_")  # by PS stratum
table(dta_ps_single$data$"_grp_",
      dta_ps_single$data$"_strata_")  # by data source and PS stratum

### Balance assessment of PS stratification.
plot(dta_ps_single, "balance")
plot(dta_ps_single, "ps")

### Obtain discounting parameters.
ps_bor_single <- rwe_ps_borrow(dta_ps_single, total_borrow = 30,
                               method = "distance", metric = "ovl")
ps_bor_single

### PSKM, single arm study, time-to-event outcome.
rst_km <- rwe_ps_survkm(ps_bor_single,
                        v_time    = "Y_Surv",
                        v_event   = "Status",
                        pred_tp  = 365)
rst_km
