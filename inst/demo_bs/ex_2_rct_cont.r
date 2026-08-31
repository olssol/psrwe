### Example of RCT and continuous outcome
suppressMessages(library(psrwe, quietly = TRUE))
org_digits <- options(digits = 3)
data(ex_dta_rct)

### Obtain PSs.
dta_ps_rct <- psrwe_est(ex_dta_rct,
                        v_covs = paste("V", 1:7, sep = ""),
                        v_grp = "Group", cur_grp_level = "current",
                        v_arm = "Arm", ctl_arm_level = "control",
                        ps_method = "logistic", nstrata = 5,
                        stra_ctl_only = FALSE)

### Obtain discounting parameters.
ps_bor_rct <- psrwe_borrow(dta_ps_rct, total_borrow = 30)

### PSCL, two-arm RCT, continuous outcome.
rst_cl_rct <- psrwe_compl(ps_bor_rct,
                          outcome_type = "continuous",
                          v_outcome = "Y_Con")

### Outcome analysis.
oa_cl_rct <- psrwe_outana(rst_cl_rct, alternative = "greater")
print(oa_cl_rct, show_rct = TRUE)

### Use simple Bootstrap stderr. This may take a while longer.
set.seed(12341)
rst_cl_rct_sbs <- psrwe_compl(ps_bor_rct,
                              outcome_type = "continuous",
                              v_outcome = "Y_Con",
                              stderr_method = "sbs")
oa_cl_rct_sbs <- psrwe_outana(rst_cl_rct_sbs, alternative = "greater")
print(oa_cl_rct_sbs, show_rct = TRUE)

### Use complex Bootstrap stderr. This may take a while longer.
set.seed(12342)
rst_cl_rct_cbs <- psrwe_compl(ps_bor_rct,
                              outcome_type = "continuous",
                              v_outcome = "Y_Con",
                              stderr_method = "cbs")
oa_cl_rct_cbs <- psrwe_outana(rst_cl_rct_cbs, alternative = "greater")
print(oa_cl_rct_cbs, show_rct = TRUE)

## Reset to user's options.
options(org_digits)

