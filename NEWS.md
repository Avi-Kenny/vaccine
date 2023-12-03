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
