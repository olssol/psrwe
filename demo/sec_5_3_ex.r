### Example of Section 5.3.
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

### PSCL, single arm study, binary outcome, weights of ATT.
rst_cl <- psrwe_compl(ps_bor_single,
                      outcome_type = "binary",
                      v_outcome    = "Y_Bin")
rst_cl


### Outcome analysis.
oa_cl <- psrwe_outana(rst_cl, mu = 0.4)
oa_cl

