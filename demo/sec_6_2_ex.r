### Example of Section 5.2.
suppressMessages(library(psrwe, quietly = TRUE))
org_digits <- options(digits = 3)
data(ex_dta_rct)
ex_dta_rct$Y_Bin <- ifelse(ex_dta_rct$Y_Con < 320, 1, 0)

### First parts of Data.
head(ex_dta_rct)

### Obtain PSs.
dta_ps_rct <- psrwe_est(ex_dta_rct,
                        v_covs = paste("V", 1:7, sep = ""),
                        v_grp = "Group", cur_grp_level = "current",
                        v_arm = "Arm", ctl_arm_level = "control",
                        ps_method = "logistic", nstrata = 1,
                        trim_ab = "none")

### Obtain discounting parameters.
ps_bor_rct <- psrwe_borrow(dta_ps_rct, total_borrow = 30)
ps_bor_xaworg <- psrwe_borrow(dta_ps_rct, total_borrow = 1)  # Original


### PSPP, RCT, binary outcome, weights of ATT.
org_mc.cores <- options(mc.cores = 1)
.msg <- capture.output({ suppressWarnings({
rst_pp_rct <- psrwe_powerp_watt(ps_bor_rct,
                                outcome_type = "binary",
                                v_outcome    = "Y_Bin",
                                seed         = 1234)
}) })
rst_pp_rct

### PSPP, RCT, binary outcome, weights of ATT, Xi.Ada.Wang.
org_mc.cores <- options(mc.cores = 1)
.msg <- capture.output({ suppressWarnings({
rst_pp_xaw <- psrwe_powerp_watt(ps_bor_rct,
                                outcome_type = "binary",
                                ipw_method   = "Xi.Ada.Wang",
                                v_outcome    = "Y_Bin",
                                seed         = 1234)
}) })
rst_pp_xaw

### PSPP, RCT, binary outcome, weights of ATT, Xi.Ada.Wang original.
org_mc.cores <- options(mc.cores = 1)
.msg <- capture.output({ suppressWarnings({
rst_pp_xaworg <- psrwe_powerp_watt(ps_bor_xaworg,                 # A = 1
                                   outcome_type = "binary",
                                   ipw_method   = "Xi.Ada.Wang",  # Original
                                   v_outcome    = "Y_Bin",
                                   seed         = 1234)
}) })
rst_pp_xaw


### Outcome analysis.
oa_pp_rct <- psrwe_outana(rst_pp_rct, alternative = "greater")
print(oa_pp_rct, show_rct = TRUE)
oa_pp_xaw <- psrwe_outana(rst_pp_xaw, alternative = "greater")
print(oa_pp_xaw, show_rct = TRUE)
oa_pp_xaworg <- psrwe_outana(rst_pp_xaworg, alternative = "greater")
print(oa_pp_xaworg, show_rct = TRUE)

### Check for A roughly in Xi.Ada.Wang original.
eps <- ps_bor_rct$data$"_ps_"[ps_bor_rct$data$"_grp_" == 0]
A_xaworg <- 1 / mean(eps / (1 - eps))
print(A_xaworg)

### Reset to user's options.
options(c(org_digits, org_mc.cores))

