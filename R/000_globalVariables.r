## Do not delete this file and file name!!
## This file should be loaded before all other *.r files.

## This is to avoid the false positive messages from R CMD check.
## "no visible binding for global variable"
## Suggested by Prof Brian Ripley
## ? globalVariables

if (getRversion() >= "2.15.1") {
    utils::globalVariables(c("_grp_",       "_arm_",   "_strata_",
                             "Stratum",     "Group",   "theta",
                             "Arm",         "Y",       "label",
                             "_matchid_",   "_ps_",    "dif_ps",
                             "_matchn_",    "asd",     "Design",
                             "Strata",      "Lower",   "Upper",
                             "Mean",        "Arm_by_Stratum"))
}
