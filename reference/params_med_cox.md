# Set parameters controlling Cox model estimation of mediation effects

This should be used in conjunction with
[`est_med`](https://avi-kenny.github.io/vaccine/reference/est_med.md) to
set parameters controlling Cox model estimation of mediation effects;
see examples.

## Usage

``` r
params_med_cox(edge_ind = F)
```

## Arguments

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
ests_med <- est_med(
  dat = dat,
  type = "Cox",
  t_0 = 578,
  params_np = params_med_cox(edge_ind=TRUE)
)
# }
```
