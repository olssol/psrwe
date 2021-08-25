### Example of Section 4.1.
suppressMessages(library(psrwe, quietly = TRUE))
data(ex_dta)

### First parts of Data.
head(ex_dta)

### Obtain PSs.
dta_ps_single <- rwe_ps_est(ex_dta,
                     v_covs = paste("V", 1:7, sep = ""),
                     v_grp = "Group", cur_grp_level = "current",
                     ps_method = "logistic")

### PS matching.
dta_ps_match <- rwe_ps_match(dta_ps_single, ratio = 2, strata_covs = "V1")
dta_ps_match

### Balance assessment of PS stratification.
plot(dta_ps_match, "balance")
plot(dta_ps_match, "ps")
plot(dta_ps_match, "diff")
plot(dta_ps_match, "diff", metric = "astd", avg_only = TRUE)

### Obtain discounting parameters.
ps_bor_match <- rwe_ps_borrow(dta_ps_match, total_borrow = 30)
ps_bor_match

### PSCL, single arm study, binary outcome.
rst_cl <- rwe_ps_compl(ps_bor_match,
                       outcome_type = "binary",
                       v_outcome    = "Y_Bin")
rst_cl

### 95% two-sided CI
rst_cl_wl <- rwe_ps_ci(rst_cl)
rst_cl_wl$CI$Control$Overall_Estimate
rst_cl_ws <- rwe_ps_ci(rst_cl, method_ci = "wilson")
rst_cl_ws$CI$Control$Overall_Estimate
rst_cl_ws_p <- rwe_ps_ci(rst_cl, method_ci = "wilson",
                         method_stderr = "plain")
rst_cl_ws_p$CI$Control$Overall_Estimate

### Use optmatch with caliper
dta_ps_match_opt <- rwe_ps_match(dta_ps_single, ratio = 2, strata_covs = "V2",
                                 mat_method = "optm", caliper = 0.5)
ps_bor_match_opt <- rwe_ps_borrow(dta_ps_match_opt, total_borrow = 30)
rst_cl_opt <- rwe_ps_compl(ps_bor_match_opt,
                           outcome_type = "binary",
                           v_outcome    = "Y_Bin")
rst_cl_opt

