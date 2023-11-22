### Example of Section 5.3.
suppressMessages(library(psrwe, quietly = TRUE))
options(digits = 3)
data(ex_dta_rct)
ex_dta_rct$Y_Bin <- ifelse(ex_dta_rct$Y_Con < 320, 1, 0)

### First parts of Data.
head(ex_dta_rct)

### Obtain PSs.
dta_ps_rct <- psrwe_est(ex_dta_rct,
                     v_covs = paste("V", 1:7, sep = ""),
                     v_grp = "Group", cur_grp_level = "current",
                     ps_method = "logistic", nstrata = 1)

### Obtain discounting parameters.
ps_bor_rct <- psrwe_borrow(dta_ps_rct, total_borrow = 30)

### PSPP, RCT, binary outcome, weights of ATT.
options(mc.cores = 1)
.msg <- capture.output({ suppressWarnings({
rst_pp_rct <- psrwe_powerp_watt(ps_bor_rct,
                            outcome_type = "binary",
                            v_outcome    = "Y_Bin",
                            seed         = 1234)
}) })
rst_pp_rct

### PSPP, single arm study, binary outcome, weights of ATT, analytic solution.
rst_ppana_rct <- psrwe_powerp_watt(ps_bor_rct,
                            outcome_type = "binary",
                            v_outcome    = "Y_Bin",
                            mcmc_method  = "analytic",
                            seed         = 1234)
rst_ppana_rct

### Outcome analysis.
oa_pp_rct <- psrwe_outana(rst_pp_rct, mu = 0.4)
oa_pp_rct
oa_ppana_rct <- psrwe_outana(rst_ppana_rct, mu = 0.4)
oa_ppana_rct

