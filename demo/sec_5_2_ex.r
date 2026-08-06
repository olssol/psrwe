### Example of Section 5.2.
suppressMessages(library(psrwe, quietly = TRUE))
org_digits <- options(digits = 3)
data(ex_dta_rct)
ex_dta_rct$Y_Bin <- ifelse(ex_dta_rct$Y_Con < 320, 1, 0)

### First parts of Data.
head(ex_dta_rct)

### Obtain PSs.
dta_ps_rct <- psrwe_est(ex_dta_rct,
                        v_covs = paste("V", 1:7, sep = ""),
                        v_grp = "Group", cur_grp_level = "current",
                        v_arm = "Arm", ctl_arm_level = "control",
                        ps_method = "logistic", nstrata = 1,
                        trim_ab = "none")

### Obtain discounting parameters.
ps_bor_rct <- psrwe_borrow(dta_ps_rct, total_borrow = 30)


### PSPP, RCT, binary outcome, weights of ATT.
org_mc.cores <- options(mc.cores = 1)
.msg <- capture.output({ suppressWarnings({
rst_pp_rct <- psrwe_powerp_watt(ps_bor_rct,
                                outcome_type = "binary",
                                v_outcome    = "Y_Bin",
                                seed         = 1234)
}) })
rst_pp_rct

### PSPP, RCT, binary outcome, weights of ATT, analytic solution.
rst_ppana_rct <- psrwe_powerp_watt(ps_bor_rct,
                                   outcome_type = "binary",
                                   v_outcome    = "Y_Bin",
                                   mcmc_method  = "analytic",
                                   seed         = 1234)
rst_ppana_rct

### Outcome analysis.
oa_pp_rct <- psrwe_outana(rst_pp_rct, alternative = "greater")
print(oa_pp_rct, show_rct = TRUE)
oa_ppana_rct <- psrwe_outana(rst_ppana_rct, alternative = "greater")
print(oa_ppana_rct, show_rct = TRUE)


### PSPP, RCT, continuous outcome, weights of ATT.
org_mc.cores <- options(mc.cores = 1)
.msg <- capture.output({ suppressWarnings({
rst_pp_rct_con <- psrwe_powerp_watt(ps_bor_rct,
                                    outcome_type = "continuous",
                                    v_outcome    = "Y_Con",
                                    seed         = 1234)
}) })
rst_pp_rct_con

### PSPP, RCT, continuous outcome, weights of ATT, analytic solution.
rst_ppana_rct_con <- psrwe_powerp_watt(ps_bor_rct,
                                       outcome_type = "continuous",
                                       v_outcome    = "Y_Con",
                                       mcmc_method  = "analytic",
                                       seed         = 1234)
rst_ppana_rct_con

### Outcome analysis.
oa_pp_rct_con <- psrwe_outana(rst_pp_rct_con, alternative = "greater")
print(oa_pp_rct_con, show_rct = TRUE)
oa_ppana_rct_con <- psrwe_outana(rst_ppana_rct_con, alternative = "greater")
print(oa_ppana_rct_con, show_rct = TRUE)


### PSPP, RCT, continuous outcome, weights of ATT, weighted for tau0.
.msg <- capture.output({ suppressWarnings({
rst_ppwtau_rct_con <- psrwe_powerp_watt(ps_bor_rct,
                                        outcome_type = "continuous",
                                        v_outcome    = "Y_Con",
                                        tau0_method  = "weighted",
                                        seed         = 1234)
}) })
rst_ppwtau_rct_con

### PSPP, RCT, continuous outcome, weights of ATT, analytic solution,
### weighted for tau0.
rst_ppanawtau_rct_con <- psrwe_powerp_watt(ps_bor_rct,
                                           outcome_type = "continuous",
                                           v_outcome    = "Y_Con",
                                           mcmc_method  = "analytic",
                                           tau0_method  = "weighted",
                                           seed         = 1234)
rst_ppanawtau_rct_con

### Outcome analysis.
oa_ppwtau_rct_con <- psrwe_outana(rst_ppwtau_rct_con, alternative = "greater")
print(oa_ppwtau_rct_con, show_rct = TRUE)
oa_ppanawtau_rct_con <- psrwe_outana(rst_ppanawtau_rct_con, alternative = "greater")
print(oa_ppanawtau_rct_con, show_rct = TRUE)


### PSPP, single arm study, continuous outcome, weights of ATT, wattcon.
.msg <- capture.output({ suppressWarnings({
rst_pp_rct_wattcon <- psrwe_powerp_watt(ps_bor_rct,
                                        outcome_type = "continuous",
                                        v_outcome    = "Y_Con",
                                        mcmc_method  = "wattcon",
                                        seed         = 1234)
}) })
rst_pp_rct_wattcon

### Outcome analysis.
oa_pp_rct_wattcon <- psrwe_outana(rst_pp_rct_wattcon, alternative = "greater")
print(oa_pp_rct_wattcon, show_rct = TRUE)

### Reset to user's options.
options(c(org_digits, org_mc.cores))

