### Example of Section 4.3.
suppressMessages(library(psrwe, quietly = TRUE))
org_options <- options(digits = 3)
data(ex_dta)

### First parts of Data.
head(ex_dta)

### Obtain PSs.
dta_ps_single <- psrwe_est(ex_dta,
                     v_covs = paste("V", 1:7, sep = ""),
                     v_grp = "Group", cur_grp_level = "current",
                     ps_method = "logistic")

### PS matching.
dta_ps_match <- psrwe_match(dta_ps_single, ratio = 2, strata_covs = "V1",
                            seed = 123)
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

### Outcome analysis.
oa_cl <- psrwe_outana(rst_cl, method_ci = "wilson", mu = 0.40)
oa_cl

### Use optmatch with caliper.
.remark_optmatch <- function() {
    warning("The optmatch may restrict use (academic license).")
    dta_ps_match_opt <- psrwe_match(dta_ps_single, ratio = 2,
                                    strata_covs = "V2",
                                    mat_method = "optm", caliper = 0.5,
                                    seed = 123)
    ps_bor_match_opt <- psrwe_borrow(dta_ps_match_opt, total_borrow = 30)
    rst_cl_opt <- psrwe_compl(ps_bor_match_opt,
                              outcome_type = "binary",
                              v_outcome    = "Y_Bin")
    rst_cl_opt
}
### Unmark below to run if optmatch is available.
# .remark_optmatch()

### Reset to user's options.
options(org_digits)

