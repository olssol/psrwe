### Example of Section 4.1.
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

### PSPP, single arm study, binary outcome, weights of ATT.
options(mc.cores = 1)
.msg <- capture.output({ suppressWarnings({
rst_pp_ns5 <- psrwe_powerp_watt(ps_bor_single_ns5,
                            outcome_type = "binary",
                            v_outcome    = "Y_Bin",
                            seed         = 1234)
}) })
rst_pp_ns5

.msg <- capture.output({ suppressWarnings({
rst_pp_ns1 <- psrwe_powerp_watt(ps_bor_single_ns1,
                            outcome_type = "binary",
                            v_outcome    = "Y_Bin",
                            seed         = 1234)
}) })
rst_pp_ns1

### PSPP, single arm study, binary outcome (from sec_4_1_ex)
.msg <- capture.output({ suppressWarnings({
rst_pp <- psrwe_powerp(ps_bor_single_ns5,
                       outcome_type = "binary",
                       v_outcome    = "Y_Bin",
                       seed         = 1234)
}) })
rst_pp

### Outcome analysis.
oa_pp_ns5 <- psrwe_outana(rst_pp_ns5, mu = 0.4)
oa_pp_ns5
oa_pp_ns1 <- psrwe_outana(rst_pp_ns1, mu = 0.4)
oa_pp_ns1
oa_pp <- psrwe_outana(rst_pp, mu = 0.4)
oa_pp

