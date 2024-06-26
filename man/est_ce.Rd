% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/est_ce.R
\name{est_ce}
\alias{est_ce}
\title{Estimate controlled effect curves}
\usage{
est_ce(
  dat,
  type = "Cox",
  t_0,
  cr = TRUE,
  cve = FALSE,
  cr_placebo_arm = F,
  s_out = seq(from = min(dat$s, na.rm = TRUE), to = max(dat$s, na.rm = TRUE), l = 101),
  ci_type = "transformed",
  placebo_risk_method = "KM",
  return_p_value = FALSE,
  return_extras = FALSE,
  params_cox = params_ce_cox(),
  params_np = params_ce_np()
)
}
\arguments{
\item{dat}{A data object returned by load_data}

\item{type}{One of c("Cox", "NP"). This specifies whether to estimate the
curve(s) using a marginalized Cox proportional hazards model or using a
monotone-constrained nonparametric estimator.}

\item{t_0}{Time point of interest}

\item{cr}{Boolean. If TRUE, the controlled risk (CR) curve is computed and
returned.}

\item{cve}{Boolean. If TRUE, the controlled vaccine efficacy (CVE) curve is
computed and returned.}

\item{cr_placebo_arm}{Boolean. If TRUE, the CR curve is estimated for the
placebo arm instead of the vaccine arm.}

\item{s_out}{A numeric vector of s-values (on the biomarker scale) for which
cve(s) and/or cr(s) are computed. Defaults to a grid of 101 points
between the min and max biomarker values.}

\item{ci_type}{One of c("transformed", "truncated", "regular", "none"). If
ci_type="transformed", confidence intervals are computed on the logit(CR)
and/or log(1-CVE) scale to ensure that confidence limits lie within [0,1]
for CR and/or lie within (-inf,1] for CVE. If ci_type="truncated",
confidence limits are constructed on the CR and/or CVE scale but
truncated to lie within [0,1] for CR and/or lie within (-inf,1] for CVE.
If ci_type="regular", confidence limits are not transformed or truncated.
If ci_type="none", confidence intervals are not computed.}

\item{placebo_risk_method}{One of c("KM", "Cox"). Method for estimating
overall risk in the placebo group. "KM" computes a Kaplan-Meier estimate
and "Cox" computes an estimate based on a marginalized Cox model survival
curve. Only relevant if cve=TRUE.}

\item{return_p_value}{Boolean; if TRUE, a P-value corresponding to the null
hypothesis that the CVE curve is flat is returned. The type of P-value
corresponds to the \code{type} argument.}

\item{return_extras}{Boolean; if TRUE, objects useful for debugging are
returned.}

\item{params_cox}{A list of options returned by
\code{\link{params_ce_cox}} that are relevant if type="Cox".}

\item{params_np}{A list of options returned by \code{\link{params_ce_np}}
that are relevant if type="NP".}
}
\value{
A list of the form \code{list(cr=list(...), cve=list(...))}
    containing CR and/or CVE estimates. Each of the inner lists contains the
    following: \itemize{
        \item{\code{s}: a vector of marker values corresponding to s_out}
        \item{\code{est}: a vector of point estimates}
        \item{\code{ci_lower}: a vector of confidence interval lower limits}
        \item{\code{ci_upper}: a vector of confidence interval upper limits}
}
}
\description{
Estimate controlled risk (CR) curves and/or controlled vaccine
    efficacy (CVE) curves. See references for definitions of these curves.
}
\examples{
data(hvtn505)
dat <- load_data(time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt",
                 marker="IgG_V2", covariates=c("age","BMI","bhvrisk"),
                 weights="wt", ph2="casecontrol", data=hvtn505)
\donttest{
ests_cox <- est_ce(dat=dat, type="Cox", t_0=578)
ests_np <- est_ce(dat=dat, type="NP", t_0=578)
}
}
\references{
Gilbert P, Fong Y, Kenny A, and Carone, M (2022). A Controlled
    Effects Approach to Assessing Immune Correlates of Protection.
    <doi:10.1093/biostatistics/kxac024>
}
