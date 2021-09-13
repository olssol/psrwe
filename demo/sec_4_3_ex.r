### Example of Section 4.3.
suppressMessages(library(psrwe, quietly = TRUE))
options(digits = 3)
data(ex_dta)

### First parts of Data.
head(ex_dta)

### Obtain PSs.
dta_ps_single <- psrwe_est(ex_dta,
                     v_covs = paste("V", 1:7, sep = ""),
                     v_grp = "Group", cur_grp_level = "current",
                     ps_method = "logistic")

### PS matching.
dta_ps_match <- psrwe_match(dta_ps_single, ratio = 2, strata_covs = "V1")
dta_ps_match

### Balance assessment of PS stratification.
plot(dta_ps_match, "balance")
plot(dta_ps_match, "ps")
plot(dta_ps_match, "diff")
plot(dta_ps_match, "diff", metric = "astd", avg_only = TRUE)

### Obtain discounting parameters.
ps_bor_match <- psrwe_borrow(dta_ps_match, total_borrow = 30)
ps_bor_match

### PSCL, single arm study, binary outcome.
rst_cl <- psrwe_compl(ps_bor_match,
                       outcome_type = "binary",
                       v_outcome    = "Y_Bin")
rst_cl

### 95% two-sided CI
rst_cl_wl <- psrwe_ci(rst_cl)
rst_cl_wl
rst_cl <- psrwe_ci(rst_cl, method_ci = "wilson")
rst_cl
rst_cl_ws_p <- psrwe_ci(rst_cl, method_ci = "wilson",
                         method_stderr = "plain")
rst_cl_ws_p

### Inference.
rst_cl <- psrwe_infer(rst_cl, mu = 0.40)
rst_cl

### Outcome analysis.
oa_cl <- psrwe_outana(rst_cl)
oa_cl

### Use optmatch with caliper
dta_ps_match_opt <- psrwe_match(dta_ps_single, ratio = 2, strata_covs = "V2",
                                 mat_method = "optm", caliper = 0.5)
ps_bor_match_opt <- psrwe_borrow(dta_ps_match_opt, total_borrow = 30)
rst_cl_opt <- psrwe_compl(ps_bor_match_opt,
                           outcome_type = "binary",
                           v_outcome    = "Y_Bin")
rst_cl_opt

