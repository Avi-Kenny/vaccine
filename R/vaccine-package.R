#' @keywords internal
"_PACKAGE"

#' HVTN 505 Dataset
#'
#' A dataset from the HVTN 505 clinical trial
#'
#' @format A data frame with 1,950 rows and 10 variables:
#' \itemize{
#'   \item{pub_id}{x}
#'   \item{trt}{Treatment Assignment: 1 = Vaccine, 0 = Placebo}
#'   \item{age}{Age (in years) at randomization}
#'   \item{BMI}{Body Mass Index:  (Weight in Kg)/(Height in meters)**2}
#'   \item{bhvrisk}{Baseline behavioral risk score}
#'   \item{HIVwk28fu}{TO DO}
#'   \item{HIVwk28preunbl}{Indicator of HIV-1 infection diagnosis on/after study week 28 (day 196) prior to Unblinding Date (22Apr2013).}
#'   \item{casecontrol}{Indicator of inclusion in the case-control cohort}
#'   \item{wt}{inverse-sampling weights, for the case-control cohort}
#'   \item{logpctpos_scaled}{TO DO}
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
