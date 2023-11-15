### Example of Section 5.3.
suppressMessages(library(psrwe, quietly = TRUE))
options(digits = 3)
data(ex_dta)

### First parts of Data.
head(ex_dta)

### Obtain PSs.
dta_ps_single_ns5 <- psrwe_est(ex_dta,
                     v_covs = paste("V", 1:7, sep = ""),
                     v_grp = "Group", cur_grp_level = "current",
                     ps_method = "logistic", nstrata = 5)
dta_ps_single_ns1 <- psrwe_est(ex_dta,
                     v_covs = paste("V", 1:7, sep = ""),
                     v_grp = "Group", cur_grp_level = "current",
                     ps_method = "logistic", nstrata = 1)

### Obtain discounting parameters.
ps_bor_single_ns5 <- psrwe_borrow(dta_ps_single_ns5, total_borrow = 30)
ps_bor_single_ns1 <- psrwe_borrow(dta_ps_single_ns1, total_borrow = 30)

### PSCL, single arm study, binary outcome, weights of ATT.
rst_cl_ns1 <- psrwe_compl(ps_bor_single_ns1,
                      outcome_type = "binary",
                      v_outcome    = "Y_Bin")
rst_cl_ns1

### PSCL, single arm study, binary outcome (from sec_4_1_ex)
rst_cl <- psrwe_compl(ps_bor_single_ns5,
                      outcome_type = "binary",
                      v_outcome    = "Y_Bin")
rst_cl


### Outcome analysis.
oa_cl_ns1 <- psrwe_outana(rst_cl_ns1, mu = 0.4)
oa_cl_ns1
oa_cl <- psrwe_outana(rst_cl, mu = 0.4)
oa_cl

