### Example of single arm and matching
suppressMessages(library(psrwe, quietly = TRUE))
options(digits = 3)
data(ex_dta)

### Obtain PSs.
dta_ps_single <- psrwe_est(ex_dta,
                     v_covs = paste("V", 1:7, sep = ""),
                     v_grp = "Group", cur_grp_level = "current",
                     ps_method = "logistic")

### PS matching.
dta_ps_match <- psrwe_match(dta_ps_single, ratio = 2, strata_covs = "V1")

### Obtain discounting parameters.
ps_bor_match <- psrwe_borrow(dta_ps_match, total_borrow = 30)

### PSCL, single arm study, binary outcome.
rst_cl <- psrwe_compl(ps_bor_match,
                      outcome_type = "binary",
                      v_outcome    = "Y_Bin")

### Outcome analysis.
oa_cl <- psrwe_outana(rst_cl, method_ci = "wilson", mu = 0.40)
oa_cl

### Use simple Bootstrap stderr. This may take very long to finish.
rst_cl_sbs <- psrwe_compl(ps_bor_match,
                          outcome_type = "binary",
                          v_outcome = "Y_Bin",
                          stderr_method = "sbs")
oa_cl_sbs <- psrwe_outana(rst_cl_sbs, mu = 0.4)
oa_cl_sbs

### Use complex Bootstrap stderr. This may take very long to finish.
rst_cl_cbs <- psrwe_compl(ps_bor_match,
                          outcome_type = "binary",
                          v_outcome = "Y_Bin",
                          stderr_method = "cbs")
oa_cl_cbs <- psrwe_outana(rst_cl_cbs, mu = 0.4)
oa_cl_cbs

