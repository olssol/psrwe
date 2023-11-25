### Example of Section 6.1.
suppressMessages(library(psrwe, quietly = TRUE))
options(digits = 3)
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


### PSPP, single arm study, continuous outcome, weights of ATT.
options(mc.cores = 1)
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

### PSPP, single arm study, continuous outcome, weights of ATT, rstan_watt.
.msg <- capture.output({ suppressWarnings({
rst_pprstan_con <- psrwe_powerp_watt_con(ps_bor_single,
                                         outcome_type = "continuous",
                                         v_outcome    = "Y_Con",
                                         mcmc_method  = "rstan_watt",
                                         seed         = 1234)
}) })
rst_pprstan_con

### Outcome analysis.
oa_pp_con <- psrwe_outana(rst_pp_con, mu = 362)
oa_pp_con
oa_ppana_con <- psrwe_outana(rst_ppana_con, mu = 362)
oa_ppana_con
oa_pprstan_con <- psrwe_outana(rst_pprstan_con, mu = 362)
oa_pprstan_con

