### Example of Section 4.1.
suppressMessages(library(psrwe, quietly = TRUE))
data(ex_dta)

### Obtain PSs.
dta_ps <- rwe_ps(ex_dta,
                 v_covs = paste("V", 1:7, sep = ""),
                 v_grp = "Group",
                 cur_grp_level = "current",
                 nstrata = 5)

### Sample size.
print(table(dta_ps$data$"_grp_"))     # by data source
print(table(dta_ps$data$"_strata_"))  # by PS stratum
print(table(dta_ps$data$"_grp_",
            dta_ps$data$"_strata_"))  # by data source and PS stratum

### Balance assessment of PS stratification.
g.balance <- plot(dta_ps, "balance")
print(g.balance)

### PS distribution.
g.ps <- plot(dta_ps, "ps")
print(g.ps)

### Obtain discounting parameters.
ps_dist <- rwe_ps_dist(dta_ps)
print(ps_dist)

### PSPP, single arm study, binary outcome.
.msg <- capture.output({ suppressWarnings({
post_smps <- rwe_ps_powerp(dta_ps,
                           total_borrow  = 40,
                           v_distance    = ps_dist$Dist[1:dta_ps$nstrata],
                           outcome_type  = "binary",
                           v_outcome     = "Y")
}) })

### PSPP results.
res_pspp <- summary(post_smps)
print(res_pspp)

### PSCL, single arm study, binary outcome.
ps_borrow <- rwe_ps_borrow(total_borrow = 40, ps_dist)
print(ps_borrow)
rst_cl <- rwe_ps_cl(dta_ps, v_borrow = ps_borrow, v_outcome = "Y")
res_pscl <- summary(rst_cl)
print(res_pscl)

