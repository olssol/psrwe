### Example of Section 5.1.
suppressMessages(library(psrwe, quietly = TRUE))
org_digits <- options(digits = 3)
data(ex_dta)

### First parts of Data.
head(ex_dta)

### Obtain PSs.
dta_ps_single <- psrwe_est(ex_dta,
                           v_covs = paste("V", 1:7, sep = ""),
                           v_grp = "Group", cur_grp_level = "current",
                           ps_method = "logistic", nstrata = 1,
                           trim_ab = "none")

### Obtain discounting parameters.
ps_bor_single <- psrwe_borrow(dta_ps_single, total_borrow = 30)
ps_bor_xaworg <- psrwe_borrow(dta_ps_single, total_borrow = 1)    # Original


### PSPP, single arm study, binary outcome, weights of ATT.
org_mc.cores <- options(mc.cores = 1)
.msg <- capture.output({ suppressWarnings({
rst_pp <- psrwe_powerp_watt(ps_bor_single,
                            outcome_type = "binary",
                            ipw_method   = "Heng.Li",
                            v_outcome    = "Y_Bin",
                            seed         = 1234)
}) })
rst_pp

### PSPP, single arm study, binary outcome, weights of ATT, Xi.Ada.Wang.
org_mc.cores <- options(mc.cores = 1)
.msg <- capture.output({ suppressWarnings({
rst_xaw <- psrwe_powerp_watt(ps_bor_single,
                             outcome_type = "binary",
                             ipw_method   = "Xi.Ada.Wang",
                             v_outcome    = "Y_Bin",
                             seed         = 1234)
}) })
rst_xaw

### PSPP, single arm study, binary outcome, weights of ATT, Xi.Ada.Wang original.
org_mc.cores <- options(mc.cores = 1)
.msg <- capture.output({ suppressWarnings({
rst_xaworg <- psrwe_powerp_watt(ps_bor_xaworg,                 # A = 1
                                outcome_type = "binary",
                                ipw_method   = "Xi.Ada.Wang",  # Original
                                v_outcome    = "Y_Bin",
                                seed         = 1234)
}) })
rst_xaw

### Outcome analysis.
oa_pp <- psrwe_outana(rst_pp, mu = 0.4)
oa_pp
oa_xaw <- psrwe_outana(rst_xaw, mu = 0.4)
oa_xaw
oa_xaworg <- psrwe_outana(rst_xaworg, mu = 0.4)
oa_xaworg

### Check for A roughly in Xi.Ada.Wang original.
eps <- ps_bor_single$data$"_ps_"[ps_bor_single$data$"_grp_" == 0]
A_xaworg <- 1 / mean(eps / (1 - eps))
print(A_xaworg)

### Reset to user's options.
options(c(org_digits, org_mc.cores))

