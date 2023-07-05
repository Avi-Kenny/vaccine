#' Plotting controlled effect curves
#'
#' @description Plot CR and/or CVE curves
#' @param ... One or more objects of class \code{"vaccine_est"} returned by
#'     \code{\link{est_ce}}.
#' @param which One of c("CR", "CVE"); controls whether to plot CR curves or CVE
#'     curves.
#' @param labels A character vector of labels corresponding to the estimate
#'     objects.
#' @param density A list with keys \code{s} and \code{weights} used to construct
#'     and plot a kernel estimate of the marker density. \code{s} gives a vector
#'     of marker values and \code{weights} gives the corresponding
#'     inverse-probability-of-sampling weights.
#' @examples
#' \dontrun{
#' data(hvtn505)
#' dat <- load_data(time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt",
#'                  marker="logpctpos_scaled", covariates=c("age","BMI","bhvrisk"),
#'                  weights="wt", ph2="casecontrol", data=hvtn505)
#' ests_cox <- est_ce(dat=dat, type="Cox", t_0=578)
#' ests_np <- est_ce(dat=dat, type="NP", t_0=578)
#' plot_ce(ests_cox, ests_np)
#' plot_ce(ests_cox, ests_np, density=list(s=dat$v$s, weights=dat$v$weights))
#' }
#' @export
plot_ce <- function(..., which, labels=NA, density=NA) {

  # To prevent R CMD CHECK notes
  x <- y <- ci_lower <- ci_upper <- NULL; rm(x,y,ci_lower,ci_upper);

  df_plot <- data.frame(
    x = double(),
    y = double(),
    ci_lower = double(),
    ci_upper = double()
  )

  if (missing(labels)) {
    labels <- sapply(substitute(list(...))[-1], deparse)
  }

  counter <- 1
  for (obj in list(...)) {
    df_add <- data.frame(
      x = obj$cr$s,
      y = obj$cr$est,
      ci_lower = obj$cr$ci_lower,
      ci_upper = obj$cr$ci_upper,
      which = labels[counter]
    )
    df_plot <- rbind(df_plot,df_add)
    counter <- counter + 1
  }

  plot <- ggplot2::ggplot(df_plot, ggplot2::aes(x=x, y=y, color=which, fill=which))

  if (!missing(density)) {
    plot <- plot +
      ggplot2::geom_density(
        ggplot2::aes(x=x),
        data = data.frame(x=density$s),
        inherit.aes = F,
        fill = "forestgreen",
        alpha = 0.3,
        color = NA
      )
    # plot
    # inds <- which(!is.na(density$s))
    # dens <- density(x=density$s[inds], weights=density$weights[inds])
    # df_plot2 <- data.frame(x=dens$x, y=dens$y, color=NA, fill=NA)
    # plot <- plot + geom_area(data=df_plot2, fill="forestgreen", alpha=0.5, color="white")
  }

  plot <- plot + ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=ci_lower, ymax=ci_upper), # color="darkorchid3", fill="darkorchid3"
                         alpha = 0.1, linetype = "dotted") +
    ggplot2::labs(title="Controlled Risk", x="log(S)", y="Risk", color="",
                  fill="") +
    ggplot2::theme(
      legend.position = "bottom"
    )
  #

  # !!!!! Cut off at quantiles

  return(plot)

}
