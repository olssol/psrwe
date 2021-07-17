### Example of Section 4.2.
suppressMessages(library(psrwe, quietly = TRUE))
data(ex_dta_rct)

### Obtain PSs.
dta_ps_2arm <- rwe_ps(ex_dta_rct,
                      v_covs = paste("V", 1:7, sep = ""),
                      v_grp = "Group",
                      cur_grp_level = "current",
                      nstrata = 5)

### Sample size.
table(dta_ps_2arm$data$"Arm")       # by arm
table(dta_ps_2arm$data$"_grp_")     # by data source
table(dta_ps_2arm$data$"_strata_")  # by PS stratum
table(dta_ps_2arm$data$"Arm",
      dta_ps_2arm$data$"_strata_",
      dta_ps_2arm$data$"_grp_")     # by arm, PS stratum, data source

### Balance assessment of PS stratification.
plot(dta_ps_2arm, "balance")
plot(dta_ps_2arm, "ps")

### PSCL, two-arm RCT, continuous outcome.
rst_2arm <- rwe_ps_cl2arm(dta_ps_2arm,
                          v_arm = "Arm",
                          trt_arm_level = 1,
                          outcome_type = "continuous",
                          v_outcome = "Y",
                          total_borrow = 40)

### Results.
rst_2arm

