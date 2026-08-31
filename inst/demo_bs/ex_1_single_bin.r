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

### Use simple Bootstrap stderr. This may take a while longer.
set.seed(12341)
rst_cl_sbs <- psrwe_compl(ps_bor_single,
                          outcome_type = "binary",
                          v_outcome = "Y_Bin",
                          stderr_method = "sbs")
oa_cl_sbs <- psrwe_outana(rst_cl_sbs, mu = 0.4)
oa_cl_sbs

### Use complex Bootstrap stderr. This may take a while longer.
set.seed(12342)
rst_cl_cbs <- psrwe_compl(ps_bor_single,
                          outcome_type = "binary",
                          v_outcome = "Y_Bin",
                          stderr_method = "cbs")
oa_cl_cbs <- psrwe_outana(rst_cl_cbs, mu = 0.4)
oa_cl_cbs

## Reset to user's options.
options(org_digits)

