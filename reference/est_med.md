# Estimate mediation effects

Estimate mediation effects, including the natural direct effect (NDE),
the natural indirect effect (NIE), and the proportion mediated (PM). See
references for definitions of these objects.

## Usage

``` r
est_med(
  dat,
  type = "NP",
  t_0,
  nde = TRUE,
  nie = TRUE,
  pm = TRUE,
  scale = "RR",
  params_cox = params_med_cox(),
  params_np = params_med_np()
)
```

## Arguments

- dat:

  A data object returned by load_data

- type:

  One of c("NP", "Cox"). This specifies whether to estimate the effects
  using a marginalized Cox proportional hazards model or using a
  nonparametric estimator.

- t_0:

  Time point of interest

- nde:

  Boolean. If TRUE, the natural direct effect is computed and returned.

- nie:

  Boolean. If TRUE, the natural indirect effect is computed and
  returned.

- pm:

  Boolean. If TRUE, the proportion mediated is computed and returned.

- scale:

  One of c("RR", "VE"). This determines whether NDE and NIE estimates
  and CIs are computed on the risk ratio (RR) scale or the vaccine
  efficacy (VE) scale. The latter equals one minus the former.

- params_cox:

  A list of options returned by
  [`params_med_cox`](https://avi-kenny.github.io/vaccine/reference/params_med_cox.md)
  that are relevant if type="Cox".

- params_np:

  A list of options returned by
  [`params_med_np`](https://avi-kenny.github.io/vaccine/reference/params_med_np.md)
  that are relevant if type="NP".

## Value

A dataframe containing the following columns:

- `effect`: one of c("NDE", "NIE", "PM")

- `est`: point estimate

- `se`: standard error of point estimate

- `ci_lower`: a confidence interval lower limit

- `ci_upper`: a confidence interval upper limit

## References

Fay MP and Follmann DA (2023). Mediation Analyses for the Effect of
Antibodies in Vaccination \<doi:10.48550/arXiv.2208.06465\>

## Examples

``` r
data(hvtn505)
dat <- load_data(time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt",
                 marker="IgG_V2", covariates=c("age","BMI","bhvrisk"),
                 weights="wt", ph2="casecontrol", data=hvtn505)
# \donttest{
ests_cox <- est_med(dat=dat, type="Cox", t_0=578)
ests_np <- est_med(dat=dat, type="NP", t_0=578)
# }
```
