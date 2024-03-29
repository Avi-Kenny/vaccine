#' Run diagnostics
#'
#' @description Run a set of diagnostic plots. Note that for this function to
#'     work, \code{\link{est_ce}} must be run with \code{return_extras=T}.
#' @param obj An object of class \code{vaccine_est} returned by
#'     \code{\link{est_ce}}
#' @return A combined plot of model diagnostics
#' @examples
#' data(hvtn505)
#' dat <- load_data(time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt",
#'                  marker="IgG_V2", covariates=c("age","BMI","bhvrisk"),
#'                  weights="wt", ph2="casecontrol", data=hvtn505)
#' \donttest{
#' ests_np <- est_ce(dat=dat, type="NP", t_0=578, return_extras=TRUE)
#' diagnostics(ests_np)
#' }
#' @export
diagnostics <- function(obj) {

  # To prevent R CMD CHECK notes
  s <- est <- ind <- NULL; rm(s,est,ind);

  if (!methods::is(obj,"vaccine_est")) {
    stop(paste0("`obj` must be an object of class 'vaccine_est' returned by es",
                "t_ce()."))
  }

  p1 <- ggplot2::ggplot(obj$extras$r_Mn, ggplot2::aes(x=s, y=est)) +
    ggplot2::geom_line(color="#56B4E9") +
    ggplot2::geom_point(color="#56B4E9", size=1) +
    ggplot2::labs(y="", x="log(S)", title="r_Mn")
  p2 <- ggplot2::ggplot(obj$extras$deriv_r_Mn, ggplot2::aes(x=s, y=est)) +
    ggplot2::geom_line(color="#56B4E9") +
    ggplot2::geom_point(color="#56B4E9", size=1) +
    ggplot2::labs(y="", x="log(S)", title="deriv_r_Mn")
  p3 <- ggplot2::ggplot(obj$extras$Gamma_os_n, ggplot2::aes(x=s, y=est)) +
    ggplot2::geom_line(color="#56B4E9") +
    ggplot2::geom_point(color="#56B4E9", size=1) +
    ggplot2::labs(y="", x="log(S)", title="Gamma_os_n")
  p4 <- ggplot2::ggplot(obj$extras$f_s_n, ggplot2::aes(x=s, y=est)) +
    ggplot2::geom_line(color="#56B4E9") +
    ggplot2::geom_point(color="#56B4E9", size=1) +
    ggplot2::labs(y="", x="log(S)", title="f_s_n")

  p5 <- ggplot2::ggplot(
    obj$extras$Q_n,
    ggplot2::aes(x=t, y=est, group=factor(interaction(ind,s)), color=factor(s))
  ) +
    ggplot2::geom_line(alpha=0.5) +
    ggplot2::facet_wrap(~s) +
    ggplot2::labs(y="", title="Q_n", color="s") +
    ggplot2::theme(legend.position="none")

  p6 <- ggplot2::ggplot(
    obj$extras$Qc_n,
    ggplot2::aes(x=t, y=est, group=factor(interaction(ind,s)), color=factor(s))
  ) +
    ggplot2::geom_line(alpha=0.5) +
    ggplot2::facet_wrap(~s) +
    ggplot2::labs(y="", title="Qc_n", color="s") +
    ggplot2::theme(legend.position="none")

  plot <- ggpubr::ggarrange(plotlist=list(p5,p1,p2,p6,p3,p4), ncol=3, nrow=2)

  return(plot)

}
