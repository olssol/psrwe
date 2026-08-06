### Example of Section 5.1.
suppressMessages(library(psrwe, quietly = TRUE))
org_digits <- options(digits = 3)
data(ex_dta)

### First parts of Data.
head(ex_dta)

### Obtain PSs.
dta_ps_single <- psrwe_est(ex_dta,
                           v_covs = paste("V", 1:7, sep = ""),
                           v_grp = "Group", cur_grp_level = "current",
                           ps_method = "logistic", nstrata = 1,
                           trim_ab = "none")

### Obtain discounting parameters.
ps_bor_single <- psrwe_borrow(dta_ps_single, total_borrow = 30)


### PSPP, single arm study, binary outcome, weights of ATT.
org_mc.cores <- options(mc.cores = 1)
.msg <- capture.output({ suppressWarnings({
rst_pp <- psrwe_powerp_watt(ps_bor_single,
                            outcome_type = "binary",
                            v_outcome    = "Y_Bin",
                            seed         = 1234)
}) })
rst_pp

### PSPP, single arm study, binary outcome, weights of ATT, analytic solution.
rst_ppana <- psrwe_powerp_watt(ps_bor_single,
                               outcome_type = "binary",
                               v_outcome    = "Y_Bin",
                               mcmc_method  = "analytic",
                               seed         = 1234)
rst_ppana

### Outcome analysis.
oa_pp <- psrwe_outana(rst_pp, mu = 0.4)
oa_pp
oa_ppana <- psrwe_outana(rst_ppana, mu = 0.4)
oa_ppana


### PSPP, single arm study, continuous outcome, weights of ATT.
org_mc.cores <- options(mc.cores = 1)
.msg <- capture.output({ suppressWarnings({
rst_pp_con <- psrwe_powerp_watt(ps_bor_single,
                                outcome_type = "continuous",
                                v_outcome    = "Y_Con",
                                seed         = 1234)
}) })
rst_pp_con

### PSPP, single arm study, continuous outcome, weights of ATT,
### analytic solution.
rst_ppana_con <- psrwe_powerp_watt(ps_bor_single,
                                   outcome_type = "continuous",
                                   v_outcome    = "Y_Con",
                                   mcmc_method  = "analytic",
                                   seed         = 1234)
rst_ppana_con

### Outcome analysis.
oa_pp_con <- psrwe_outana(rst_pp_con, mu = 362)
oa_pp_con
oa_ppana_con <- psrwe_outana(rst_ppana_con, mu = 362)
oa_ppana_con


### PSPP, single arm study, continuous outcome, weights of ATT.
.msg <- capture.output({ suppressWarnings({
rst_ppwtau_con <- psrwe_powerp_watt(ps_bor_single,
                                    outcome_type = "continuous",
                                    v_outcome    = "Y_Con",
                                    tau0_method  = "weighted",
                                    seed         = 1234)
}) })
rst_ppwtau_con

### PSPP, single arm study, continuous outcome, weights of ATT,
### analytic solution, weighted for tau0.
rst_ppanawtau_con <- psrwe_powerp_watt(ps_bor_single,
                                       outcome_type = "continuous",
                                       v_outcome    = "Y_Con",
                                       mcmc_method  = "analytic",
                                       tau0_method  = "weighted",
                                       seed         = 1234)
rst_ppanawtau_con

### Outcome analysis.
oa_ppwtau_con <- psrwe_outana(rst_ppwtau_con, mu = 362)
oa_ppwtau_con
oa_ppanawtau_con <- psrwe_outana(rst_ppanawtau_con, mu = 362)
oa_ppanawtau_con


### PSPP, single arm study, continuous outcome, weights of ATT, wattcon.
.msg <- capture.output({ suppressWarnings({
rst_pp_wattcon <- psrwe_powerp_watt(ps_bor_single,
                                    outcome_type = "continuous",
                                    v_outcome    = "Y_Con",
                                    mcmc_method  = "wattcon",
                                    seed         = 1234)
}) })
rst_pp_wattcon

### Outcome analysis.
oa_pp_wattcon <- psrwe_outana(rst_pp_wattcon, mu = 362)
oa_pp_wattcon

### Reset to user's options.
options(c(org_digits, org_mc.cores))

