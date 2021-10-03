### Example of Section 4.1.
suppressMessages(library(psrwe, quietly = TRUE))
options(digits = 3)
data(ex_dta)

### First parts of Data.
head(ex_dta)

### Obtain PSs.
dta_ps_single <- psrwe_est(ex_dta,
                     v_covs = paste("V", 1:7, sep = ""),
                     v_grp = "Group", cur_grp_level = "current",
                     ps_method = "logistic", nstrata = 5)

### Balance assessment of PS stratification.
plot(dta_ps_single, "balance")
plot(dta_ps_single, "ps")
plot(dta_ps_single, "diff")
plot(dta_ps_single, "diff", metric = "astd", avg_only = TRUE)

### Obtain discounting parameters.
ps_bor_single <- psrwe_borrow(dta_ps_single, total_borrow = 30)
ps_bor_single

### PSPP, single arm study, binary outcome.
options(mc.cores = 1)
.msg <- capture.output({ suppressWarnings({
rst_pp <- psrwe_powerp(ps_bor_single,
                       outcome_type = "binary",
                       v_outcome    = "Y_Bin",
                       seed         = 1234)
}) })
rst_pp

### Plot PSPP results.
plot(rst_pp)
plot(rst_pp, add_stratum = TRUE)

### Outcome analysis.
oa_pp <- psrwe_outana(rst_pp, mu = 0.4)
oa_pp

### PSCL, single arm study, binary outcome.
rst_cl <- psrwe_compl(ps_bor_single,
                      outcome_type = "binary",
                      v_outcome    = "Y_Bin")
rst_cl

### Outcome analysis.
oa_cl <- psrwe_outana(rst_cl, mu = 0.4)
oa_cl
