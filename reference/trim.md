# Trim data for plotting/reporting

Removes a subset of estimates returned by
[`est_ce`](https://avi-kenny.github.io/vaccine/reference/est_ce.md)

## Usage

``` r
trim(ests, dat, quantiles, placebo = FALSE)
```

## Arguments

- ests:

  An object of class `"vaccine_est"` returned by
  [`est_ce`](https://avi-kenny.github.io/vaccine/reference/est_ce.md).

- dat:

  The data object originally passed into
  [`est_ce`](https://avi-kenny.github.io/vaccine/reference/est_ce.md).

- quantiles:

  A vector of length 2 representing the quantiles of the marker
  distribution at which to trim the data; if, for example,
  `quantiles=c(0.1,0.9)` is specified, values outside the 10 (weighted)
  quantiles of the marker distribution will be trimmed.

- placebo:

  Boolean; if TRUE, quantiles are computed based on the marker
  distribution in the placebo arm instead of the vaccine arm

## Value

A modified copy of `ests` with the data trimmed.

## Examples

``` r
# \donttest{
data(hvtn505)
dat <- load_data(time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt",
                 marker="IgG_V2", covariates=c("age","BMI","bhvrisk"),
                 weights="wt", ph2="casecontrol", data=hvtn505)
ests_cox <- est_ce(dat=dat, type="Cox", t_0=578)
ests_cox_tr <- trim(ests_cox, dat=dat, quantiles=c(0.05,0.95))
plot_ce(ests_cox_tr, density_type="kde", dat=dat)
# }
```
