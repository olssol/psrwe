### Example of Section 4.2.
suppressMessages(library(psrwe, quietly = TRUE))
options(digits = 3)
data(ex_dta_rct)

### First parts of Data.
head(ex_dta_rct)

### Obtain PSs.
dta_ps_rct <- psrwe_est(ex_dta_rct,
                        v_covs = paste("V", 1:7, sep = ""),
                        v_grp = "Group", cur_grp_level = "current",
                        v_arm = "Arm", ctl_arm_level = "control",
                        ps_method = "logistic", nstrata = 5,
                        stra_ctl_only = FALSE)
dta_ps_rct

### Balance assessment of PS stratification.
plot(dta_ps_rct, "balance")
plot(dta_ps_rct, "ps")
plot(dta_ps_rct, "diff")
plot(dta_ps_rct, "diff", metric = "astd", avg_only = TRUE)

### Obtain discounting parameters.
ps_bor_rct <- psrwe_borrow(dta_ps_rct, total_borrow = 30)
ps_bor_rct

### PSCL, two-arm RCT, continuous outcome.
rst_cl_rct <- psrwe_compl(ps_bor_rct,
                          outcome_type = "continuous",
                          v_outcome = "Y_Con")
rst_cl_rct

### Outcome analysis.
oa_cl_rct <- psrwe_outana(rst_cl_rct, alternative = "greater")
print(oa_cl_rct, show_rct = TRUE)

### PSPP, two-arm RCT, continuous outcome.
options(mc.cores = 1)
.msg <- capture.output({ suppressWarnings({
rst_pp_rct <- psrwe_powerp(ps_bor_rct,
                           outcome_type = "continuous",
                           v_outcome    = "Y_Con",
                           seed         = 1234)
}) })
plot(rst_pp_rct)
plot(rst_pp_rct, split_rct_arm = TRUE)
plot(rst_pp_rct, add_stratum = TRUE)
plot(rst_pp_rct, add_stratum = TRUE, split_rct_arm = TRUE)

### Outcome analysis.
oa_pp_rct <- psrwe_outana(rst_pp_rct, alternative = "greater")
print(oa_pp_rct, show_rct = TRUE)

