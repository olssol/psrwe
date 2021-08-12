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

### Obtain discounting parameters.
# ps_bor_match_ovl <- rwe_ps_borrow(dta_ps_match, total_borrow = 30,
#                                   metric = "ovl")
# ps_bor_match_ovl
ps_bor_match <- rwe_ps_borrow(dta_ps_match, total_borrow = 30)
ps_bor_match

### PSCL, single arm study, binary outcome.
rst_cl <- rwe_ps_compl(ps_bor_match,
                       outcome_type = "binary",
                       v_outcome    = "Y_Bin")
rst_cl

