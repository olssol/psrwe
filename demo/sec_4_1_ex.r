### Example of Section 4.1.
suppressMessages(library(psrwe, quietly = TRUE))
data(ex_dta)

### First parts of Data.
head(ex_dta)

### Obtain PSs.
dta_ps_single <- rwe_ps_est(ex_dta,
                     v_covs = paste("V", 1:7, sep = ""),
                     v_grp = "Group", cur_grp_level = "current",
                     ps_method = "logistic", nstrata = 5)

### Balance assessment of PS stratification.
plot(dta_ps_single, "balance")
plot(dta_ps_single, "ps")
plot(dta_ps_single, "diff")
plot(dta_ps_single, "diff", metric = "astd")

### Obtain discounting parameters.
ps_bor_single <- rwe_ps_borrow(dta_ps_single, total_borrow = 30,
                               method = "distance", metric = "ovl")
ps_bor_single

### PSPP, single arm study, binary outcome.
options(mc.cores = 1)
.msg <- capture.output({ suppressWarnings({
rst_pp <- rwe_ps_powerp(ps_bor_single,
                        outcome_type = "binary",
                        v_outcome    = "Y_Bin",
                        seed         = 1234)
}) })
rst_pp

### PSCL, single arm study, binary outcome.
rst_cl <- rwe_ps_compl(ps_bor_single,
                       outcome_type = "binary",
                       v_outcome    = "Y_Bin")
rst_cl
