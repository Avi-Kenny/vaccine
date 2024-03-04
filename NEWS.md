# vaccine 1.2.1

### Minor changes

- Fixed a bug that displayed an incorrect package version number at startup.
- Tweaked the SuperLearner library involved in standard error estimation for `est_ce(..., type="NP")`, resulting in improved performance.
- Fixed a bug related to `plot_ce(..., density_type="kde edge")`.
- Fixed a minor bug related to the `SuperLearner` package introduced in version 1.2.0.

# vaccine 1.2.0

### Major changes

- Added the `trim` function, which can be used in conjunction with the `plot_ce` function to truncate the display of estimate objects produced by `est_ce`. Specifically, estimates with X-coordinates that lie outside specified quantiles of the observed distribution of the marker are removed. See package vignette and examples.
- Added the `as_table` function, which formats estimate objects produced using `est_ce` as a table.

### Minor changes

- Implemented various changes to the `plot_ce` function: added background kernel density plotting, changed default plot styling, fixed a bug that caused CR plots to be displayed instead of CVE plots when `plot_ce(..., which="CVE")` was called, and added to the plotting section of the main package vignette.

# vaccine 1.1.0

### Major changes

- In `est_ce`, added the option `return_p_value`; if set to `TRUE`, a P-value will be returned corresponding to the null hypothesis that the CVE curve (or equivalently, the CR curve in the vaccine group) is constant. If `type='Cox'` is specified, this will be a Wald-type test using the estimated Cox model parameter vector. If `type='NP'`, this will be a nonparametric test assuming monotonicity of the curve.

### Minor changes

- Fixed importing of factor/character variables in `load_data`.

# vaccine 1.0.0

### Major changes

- Added functions for nonparametric estimation and inference for mediation analysis parameters, including the natural direct effect, the natural indirect effect, and the proportion mediated.
- Adapted nonparametric and Cox-based CVE estimation and inference to handle settings in which some or all covariates are measured only in phase two (e.g., applicable to settings in which researchers want to control for a baseline biomarker measurement).
- Added confidence band monotonization for nonparametric CVE inference.

### Minor changes

- Various minor bug fixes and code speed improvements.
- Completed unit testing framework and increased code coverage to 67%.

# vaccine 0.1.0

### Major changes

- Initial package release

### Minor changes

- None
