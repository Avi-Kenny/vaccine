# Create table of estimates

Format estimates returned by
[`est_ce`](https://avi-kenny.github.io/vaccine/reference/est_ce.md) as a
table

## Usage

``` r
as_table(..., which = "CR", labels = NA)
```

## Arguments

- ...:

  One or more objects of class `"vaccine_est"` returned by
  [`est_ce`](https://avi-kenny.github.io/vaccine/reference/est_ce.md).

- which:

  One of c("CR", "CVE"); controls whether the table contains CR or CVE
  values.

- labels:

  A character vector of length equal to length(list(...)) representing
  curve labels

## Value

A table of CR or CVE values

## Examples

``` r
data(hvtn505)
dat <- load_data(time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt",
                 marker="IgG_V2", covariates=c("age","BMI","bhvrisk"),
                 weights="wt", ph2="casecontrol", data=hvtn505)
# \donttest{
ests_cox <- est_ce(dat=dat, type="Cox", t_0=578)
ests_np <- est_ce(dat=dat, type="NP", t_0=578)
ests_table <- as_table(ests_cox, ests_np)
head(ests_table)
# }
```
