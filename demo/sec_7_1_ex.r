### Example of Section 7.1.
suppressMessages(library(psrwe, quietly = TRUE))
org_digits <- options(digits = 3)
data(ex_dta)

### First parts of Data.
head(ex_dta)

### Obtain PSs.
dta_ps_single <- psrwe_est(ex_dta,
                     v_covs = paste("V", 1:7, sep = ""),
                     v_grp = "Group", cur_grp_level = "current",
                     ps_method = "logistic", nstrata = 5)

### Obtain discounting parameters.
ps_bor_single <- psrwe_borrow(dta_ps_single, total_borrow = 30)

### PSPP, single arm study, binary outcome, prior only.
org_mc.cores <- options(mc.cores = 1)
.msg <- capture.output({ suppressWarnings({
rst_pp <- psrwe_powerp(ps_bor_single,
                       outcome_type = "binary",
                       v_outcome    = "Y_Bin",
                       prioronly    = TRUE,      
                       seed         = 1234)
}) })
rst_pp

### Plot PSPP results.
plot(rst_pp)
plot(rst_pp, add_stratum = TRUE)

options(org_digits)

