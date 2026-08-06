### Example of Section 7.2.
suppressMessages(library(psrwe, quietly = TRUE))
org_digits <- options(digits = 3)
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

### Obtain discounting parameters.
ps_bor_rct <- psrwe_borrow(dta_ps_rct, total_borrow = 30)

### PSPP, two-arm RCT, continuous outcome, prior only.
org_mc.cores <- options(mc.cores = 1)
.msg <- capture.output({ suppressWarnings({
rst_pp_rct <- psrwe_powerp(ps_bor_rct,
                           outcome_type = "continuous",
                           v_outcome    = "Y_Con",
                           prioronly    = TRUE,
                           seed         = 1234)
}) })
plot(rst_pp_rct)
plot(rst_pp_rct, split_rct_arm = TRUE)
plot(rst_pp_rct, add_stratum = TRUE)
plot(rst_pp_rct, add_stratum = TRUE, split_rct_arm = TRUE)

### Reset to user's options.
options(c(org_digits, org_mc.cores))

