# HVTN 505 Dataset

A dataset from the HVTN 505 clinical trial. This is a filtered version
of the original dataset (which contains 2,502 rows) based on the week28
variable, which defines eligibility for inclusion in a correlates
analysis. Additionally, two rows were removed due to missing data.

## Usage

``` r
data(hvtn505)
```

## Format

A data frame with 2,302 rows and 10 variables:

- pub_id: Unique individual identifier

- trt: Treatment Assignment: 1=Vaccine, 0=Placebo

- HIVwk28preunbl: Indicator of HIV-1 infection diagnosis on/after study
  week 28 (day 196) prior to Unblinding Date (Apr 22, 2013).

- HIVwk28preunblfu: Follow-up time (in days) for HIV-1 infection
  diagnosis endpoint as of the Unblinding Date (22Apr2013) occuring
  on/after study week 28 (day 196).

- age: Age (in years) at randomization

- BMI: Body Mass Index: (Weight in Kg)/(Height in meters)\*\*2

- bhvrisk: Baseline behavioral risk score

- casecontrol: Indicator of inclusion in the case-control cohort

- wt: Inverse-probability-of-sampling weights, for the case-control
  cohort

- IgG_env: IgG Binding to gp120/140

- IgG_V2: IgG Binding to V1V2

- IgG_V3: IgG Binding to V3

## Source

<https://atlas.scharp.org/project/HVTN%20Public%20Data/HVTN%20505/begin.view>
