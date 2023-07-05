#' Plotting controlled effect curves
#'
#' @description TO DO
#' @param ... One or more objects of class \code{"vaccine_est"} returned by
#'     \code{\link{est_ce}}.
#' @param which One of c("CR", "CVE"); controls whether to plot CR curves or CVE
#'     curves.
#' @examples
#' print("to do")
#' @export
plot_ce <- function(..., which) {

  ggplot(
    data.frame(x = ests_cox$cr$s,
               y = ests_cox$cr$est,
               ci_lower = ests_cox$cr$ci_lower,
               ci_upper = ests_cox$cr$ci_upper),
    aes(x=x, y=y)
  ) +
    geom_abline(slope=0, intercept=ov[ov$group=="vaccine","est"], alpha=0.5) +
    geom_line(color="darkorchid3") +
    geom_density(aes(x=x), data=data.frame(x=dat$v$s), inherit.aes=F,
                 fill="forestgreen", color=NA, alpha=0.2) +
    geom_ribbon(aes(ymin=ci_lower, ymax=ci_upper), color="darkorchid3",
                fill="darkorchid3", alpha = 0.1, linetype = "dotted") +
    labs(title="Controlled Risk")

}
