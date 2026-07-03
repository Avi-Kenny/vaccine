# Calculate summary statistics

TO DO

## Usage

``` r
summary_stats(dat, quietly = FALSE)
```

## Arguments

- dat:

  A data object returned by \`load_data\`.

- quietly:

  Boolean. If true, output will not be printed.

## Value

A list containing values of various summary statistics.

## Examples

``` r
data(hvtn505)
dat <- load_data(time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt",
                 marker="IgG_V2", covariates=c("age","BMI","bhvrisk"),
                 weights="wt", ph2="casecontrol", data=hvtn505)
summary_stats(dat=dat)
```
