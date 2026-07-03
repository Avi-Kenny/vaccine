# Load and format data object

This function takes in user-supplied data and returns a data object that
can be read in by
[`summary_stats`](https://avi-kenny.github.io/vaccine/reference/summary_stats.md),
[`est_ce`](https://avi-kenny.github.io/vaccine/reference/est_ce.md),
[`est_med`](https://avi-kenny.github.io/vaccine/reference/est_med.md),
and other estimation functions. Data is expected to come from a vaccine
clinical trial, possibly involving two-phase sampling and possibly
including a biomarker of interest.

## Usage

``` r
load_data(
  time,
  event,
  vacc,
  marker,
  covariates,
  weights,
  ph2,
  strata = NA,
  data,
  covariates_ph2 = FALSE
)
```

## Arguments

- time:

  A character string; the name of the numeric variable representing
  observed event or censoring times.

- event:

  A character string; the name of the binary variable corresponding to
  whether the observed time represents an event time (1) or a censoring
  time (0). Either integer (0/1) or Boolean (T/F) values are allowed.

- vacc:

  A character string; the name of the binary variable denoting whether
  the individual is in the vaccine group (1) or the placebo group (0).
  Accepts either integer (0/1) or Boolean (T/F) values.

- marker:

  A character string; the name of the numeric variable of biomarker
  values.

- covariates:

  A character vector; the names of the covariate columns. Columns values
  should be either numeric, binary, or factors. Character columns will
  be converted into factors.

- weights:

  A character string; the name of the numeric variable containing
  inverse-probability-of-sampling (IPS) weights.

- ph2:

  A character string; the name of the binary variable representing
  whether the individual is in the phase-two cohort (1) or not (0).
  Accepts either integer (0/1) or Boolean (T/F) values.

- strata:

  A character string; the name of the variable containing strata
  identifiers (for two-phase sampling strata).

- data:

  A dataframe containing the vaccine trial data.

- covariates_ph2:

  A boolean; if at least one of the covariates is measured only in the
  phase-two cohort, set this to TRUE.

## Value

An object of class `vaccine_dat`.

## Examples

``` r
data(hvtn505)
dat <- load_data(time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt",
                 marker="IgG_V2", covariates=c("age","BMI","bhvrisk"),
                 weights="wt", ph2="casecontrol", data=hvtn505)
```
