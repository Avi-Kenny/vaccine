# Set parameters controlling Cox model estimation of controlled effect curves

This should be used in conjunction with
[`est_ce`](https://avi-kenny.github.io/vaccine/reference/est_ce.md) to
set parameters controlling Cox model estimation of controlled effect
curves; see examples.

## Usage

``` r
params_ce_cox(spline_df = NA, spline_knots = NA, edge_ind = FALSE)
```

## Arguments

- spline_df:

  An integer; if the marker is modeled flexibly within the Cox model
  linear predictor as a natural cubic spline, this option controls the
  degrees of freedom in the spline; knots are chosen to be equally
  spaced across the range of the marker.

- spline_knots:

  A numeric vector; as an alternative to specifying `spline_df`, the
  exact locations of the knots in the spline (including boundary knots)
  can be specified with this option.

- edge_ind:

  Boolean. If TRUE, an indicator variable corresponding to the lower
  limit of the marker will be included in the Cox model linear
  predictor.

## Value

A list of options.

## Examples

``` r
data(hvtn505)
dat <- load_data(time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt",
                 marker="IgG_V2", covariates=c("age","BMI","bhvrisk"),
                 weights="wt", ph2="casecontrol", data=hvtn505)
# \donttest{
ests_cox <- est_ce(
  dat = dat,
  type = "Cox",
  t_0 = 578,
  params_cox = params_ce_cox(spline_df=4)
)
# }
```
