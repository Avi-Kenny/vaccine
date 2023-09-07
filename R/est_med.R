#' Estimate mediation effects
#'
#' @description Estimate mediation effects, including the natural direct effect
#'     (NDE), the natural indirect effect (NIE), and the proportion mediated
#'     (PM). See references for definitions of these objects.
#' @param dat A data object returned by load_data
#' @param type One of c("NP", "Cox"). This specifies whether to estimate the
#'     effects using a marginalized Cox proportional hazards model or using a
#'     nonparametric estimator.
#' @param t_0 Time point of interest
#' @param nde Boolean. If TRUE, the natural direct effect is computed and
#'     returned.
#' @param nie Boolean. If TRUE, the natural indirect effect is computed and
#'     returned.
#' @param pm Boolean. If TRUE, the proportion mediated is computed and returned.
#' @param scale One of c("RR", "VE"). This determines whether estimates and CIs
#'     are computed on the risk ratio (RR) scale or the vaccine efficacy (VE)
#'     scale. The latter equals one minus the former.
#' @return A dataframe containing the following columns: \itemize{
#'     \item{\code{effect}: one of c("NDE", "NIE", "PM")}
#'     \item{\code{est}: point estimate}
#'     \item{\code{ci_lower}: a confidence interval lower limit}
#'     \item{\code{ci_upper}: a confidence interval upper limit}
#' }
#' @examples
#' data(hvtn505)
#' dat <- load_data(time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt",
#'                  marker="IgG_V2", covariates=c("age","BMI","bhvrisk"),
#'                  weights="wt", ph2="casecontrol", data=hvtn505)
#' \donttest{
#' ests_np <- est_med(dat=dat, type="NP", t_0=578)
#' }
#' @references Fay MP and Follmann DA (2023). Mediation Analyses for the Effect
#'     of Antibodies in Vaccination <doi:10.48550/arXiv.2208.06465>
#' @export
est_med <- function(
    dat, type="NP", t_0, nde=TRUE, nie=TRUE, pm=TRUE
    # ci_type="transformed", return_extras=FALSE,
    # params_cox=params_med_cox(), params_np=params_med_cox()
) {

  # !!!!! Need to refactor with est_ce functions

  # !!!!! Integrate params

  if (type=="Cox") { stop("Coming soon!") }

  if (type=="NP") {

    # !!!!! Check for vaccine group data for NIE
    # !!!!! Check for vaccine+placebo group data for NDE+PM

  }

  return(ests)

}
