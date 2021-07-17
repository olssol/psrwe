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
table(dta_ps$data$"_grp_")     # by data source
table(dta_ps$data$"_strata_")  # by PS stratum
table(dta_ps$data$"_grp_",
      dta_ps$data$"_strata_")  # by data source and PS stratum

### Balance assessment of PS stratification.
plot(dta_ps, "balance")
plot(dta_ps, "ps")

### Obtain discounting parameters.
ps_dist <- rwe_ps_dist(dta_ps)
ps_dist

### PSPP, single arm study, binary outcome.
.msg <- capture.output({ suppressWarnings({
post_smps <- rwe_ps_powerp(dta_ps,
                           total_borrow  = 40,
                           v_distance    = ps_dist$Dist[1:dta_ps$nstrata],
                           outcome_type  = "binary",
                           v_outcome     = "Y")
}) })

### PSPP results.
summary(post_smps)

### PSCL, single arm study, binary outcome.
ps_borrow <- rwe_ps_borrow(total_borrow = 40, ps_dist)
ps_borrow
rst_cl <- rwe_ps_cl(dta_ps, v_borrow = ps_borrow, v_outcome = "Y")
summary(rst_cl)

