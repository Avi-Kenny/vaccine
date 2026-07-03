# Run diagnostics

Run a set of diagnostic plots. Note that for this function to work,
[`est_ce`](https://avi-kenny.github.io/vaccine/reference/est_ce.md) must
be run with `return_extras=T`.

## Usage

``` r
diagnostics(obj)
```

## Arguments

- obj:

  An object of class `vaccine_est` returned by
  [`est_ce`](https://avi-kenny.github.io/vaccine/reference/est_ce.md)

## Value

A combined plot of model diagnostics

## Examples

``` r
data(hvtn505)
dat <- load_data(time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt",
                 marker="IgG_V2", covariates=c("age","BMI","bhvrisk"),
                 weights="wt", ph2="casecontrol", data=hvtn505)
# \donttest{
ests_np <- est_ce(dat=dat, type="NP", t_0=578, return_extras=TRUE)
diagnostics(ests_np)
# }
```
