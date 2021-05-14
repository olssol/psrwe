# psrwe 1.3.9999

## Bug fixes

* Fixed `metric_ovl` (called in `get_distance()` using option `metric = "ovl"`)
to ensure integration is over the proper support.

## Minor changes

* Added a `NEWS.md` file to track changes to the package.
* Minor updates to `get_distance.R` (renamed from `psrwe_balance.R`).
* Added a README file.
* Switched to `Authors@R` in DESCRIPTION file.
* Updated DESCRIPTION to point to repo.
* Delete junk files from GitHub repo.
* No longer export `tkExpRst`, `tkMakeLocal`, and `tkCallFun`, nor create .Rd
documentation files
* Fixed a minor bug in DESCRIPTION Suggests
* Minor code formatting updates
* Minor documentation updates
* Fixed a minor bug that occured if a tibble was passed to `rwe_ps()` instead of
a data frame.