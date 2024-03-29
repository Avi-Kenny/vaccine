#' @keywords internal
"_PACKAGE"

#' HVTN 505 Dataset
#'
#' A dataset from the HVTN 505 clinical trial
#'
#' @format A data frame with 1,950 rows and 10 variables:
#' \itemize{
#'   \item pub_id: Unique individual identifier
#'   \item trt: Treatment Assignment: 1=Vaccine, 0=Placebo
#'   \item HIVwk28preunbl: Indicator of HIV-1 infection diagnosis on/after study
#'       week 28 (day 196) prior to Unblinding Date (Apr 22, 2013).
#'   \item HIVwk28preunblfu: Follow-up time (in days) for HIV-1 infection
#'       diagnosis endpoint as of the Unblinding Date (22Apr2013) occuring
#'       on/after study week 28 (day 196).
#'   \item age: Age (in years) at randomization
#'   \item BMI: Body Mass Index: (Weight in Kg)/(Height in meters)**2
#'   \item bhvrisk: Baseline behavioral risk score
#'   \item casecontrol: Indicator of inclusion in the case-control cohort
#'   \item wt: Inverse-probability-of-sampling weights, for the case-control
#'       cohort
#'   \item IgG_env: IgG Binding to gp120/140
#'   \item IgG_V2: IgG Binding to V1V2
#'   \item IgG_V3: IgG Binding to V3
#' }
#' @source \url{https://atlas.scharp.org/cpas/project/HVTN\%20Public\%20Data/HVTN\%20505/begin.view}
#' @docType data
#' @keywords datasets
#' @name hvtn505
#' @usage data(hvtn505)
"hvtn505"

## usethis namespace: start
#' @importFrom SuperLearner SuperLearner
## usethis namespace: end
NULL
