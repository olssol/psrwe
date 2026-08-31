### Example of single arm and binary outcome
suppressMessages(library(psrwe, quietly = TRUE))
org_digits <- options(digits = 3)
data(ex_dta)

### Obtain PSs.
dta_ps_single <- psrwe_est(ex_dta,
                     v_covs = paste("V", 1:7, sep = ""),
                     v_grp = "Group", cur_grp_level = "current",
                     ps_method = "logistic", nstrata = 5)

### Obtain discounting parameters.
ps_bor_single <- psrwe_borrow(dta_ps_single, total_borrow = 30)

### PSCL, single arm study, binary outcome.
rst_cl <- psrwe_compl(ps_bor_single,
                      outcome_type = "binary",
                      v_outcome    = "Y_Bin")

### Outcome analysis.
oa_cl <- psrwe_outana(rst_cl, mu = 0.4)
oa_cl

### Use complex Jackknife stderr. This may take a while longer.
rst_cl_cjk <- psrwe_compl(ps_bor_single,
                          outcome_type = "binary",
                          v_outcome = "Y_Bin",
                          stderr_method = "cjk")
oa_cl_cjk <- psrwe_outana(rst_cl_cjk, mu = 0.4)
oa_cl_cjk

## Reset to user's options.
options(org_digits)

