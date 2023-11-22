### Example of Section 5.1.
suppressMessages(library(psrwe, quietly = TRUE))
options(digits = 3)
data(ex_dta)

### First parts of Data.
head(ex_dta)

### Obtain PSs.
dta_ps_single <- psrwe_est(ex_dta,
                     v_covs = paste("V", 1:7, sep = ""),
                     v_grp = "Group", cur_grp_level = "current",
                     ps_method = "logistic", nstrata = 1)

### Obtain discounting parameters.
ps_bor_single <- psrwe_borrow(dta_ps_single, total_borrow = 30)

### PSPP, single arm study, binary outcome, weights of ATT.
options(mc.cores = 1)
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

