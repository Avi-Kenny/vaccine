# Estimate overall risk and vaccine efficacy

Estimate overall risk and vaccine efficacy.

## Usage

``` r
est_overall(dat, t_0, method = "Cox", risk = TRUE, ve = TRUE)
```

## Arguments

- dat:

  A data object returned by load_data

- t_0:

  Time point of interest

- method:

  One of c("KM", "Cox"), corresponding to either a Kaplan-Meier
  estimator ("KM") or a marginalized Cox proportional hazards model
  ("Cox").

- risk:

  Boolean. If TRUE, the controlled risk (CR) curve is computed.

- ve:

  Boolean. If TRUE, the controlled vaccine efficacy (CVE) curve is
  computed.

## Value

A dataframe containing estimates

## Examples

``` r
data(hvtn505)
dat <- load_data(time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt",
                 marker="IgG_V2", covariates=c("age","BMI","bhvrisk"),
                 weights="wt", ph2="casecontrol", data=hvtn505)
est_overall(dat=dat, t_0=578, method="KM")
```
