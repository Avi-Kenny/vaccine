# Set parameters controlling nonparametric estimation of mediation effects

This should be used in conjunction with
[`est_med`](https://avi-kenny.github.io/vaccine/reference/est_med.md) to
set parameters controlling nonparametric estimation of mediation
effects; see examples.

## Usage

``` r
params_med_np(
  grid_size = list(y = 101, s = 101, x = 5),
  surv_type = "survML-G",
  density_type = "binning",
  density_bins = 15
)
```

## Arguments

- grid_size:

  A list with keys `y`, `s`, and `x`; controls the rounding of data
  values. Decreasing the grid size values results in shorter computation
  times, and increasing the values results in more precise estimates. If
  grid_size\$s=101, this means that a grid of 101 equally-spaced points
  (defining 100 intervals) will be created from min(S) to max(S), and
  each S value will be rounded to the nearest grid point. For
  grid_size\$y, a grid will be created from 0 to t_0, and then extended
  to max(Y). For grid_size\$x, a separate grid is created for each
  covariate column (binary/categorical covariates are ignored).

- surv_type:

  One of c("Cox", "survSL", "survML-G", "survML-L"); controls the method
  to use to estimate the conditional survival and conditional censoring
  functions. If type="Cox", a survival function based on a Cox
  proportional hazard model will be used. If type="survSL", the Super
  Learner method of Westling 2023 is used. If type="survML-G", the
  global survival stacking method of Wolock 2022 is used. If
  type="survML-L", the local survival stacking method of Polley 2011 is
  used.

- density_type:

  One of c("binning", "parametric"); controls the method to use to
  estimate the density ratio f(S\|X)/f(S).

- density_bins:

  An integer; if density_type="binning", the number of bins to use. If
  density_bins=0, the number of bins will be selected via
  cross-validation.

## Value

A list of options.

## Examples

``` r
data(hvtn505)
dat <- load_data(time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt",
                 marker="IgG_V2", covariates=c("age","BMI","bhvrisk"),
                 weights="wt", ph2="casecontrol", data=hvtn505)
# \donttest{
ests_med <- est_med(
  dat = dat,
  type = "NP",
  t_0 = 578,
  params_np = params_med_np(surv_type="survML-L")
)
# }
```
