#' Plotting controlled effect curves
#'
#' @description Plot CR and/or CVE curves
#' @param ... One or more objects of class \code{"vaccine_est"} returned by
#'     \code{\link{est_ce}}.
#' @param which One of c("CR", "CVE"); controls whether to plot CR curves or CVE
#'     curves.
#' @param labels A character vector of labels corresponding to the estimate
#'     objects.
#' @param density_type One of c("none", "kde", "kde edge", "histogram").
#'     Controls the type of estimator used for the background marker density
#'     plot. For "none", no density plot is displayed. For "kde", a weighted
#'     kernel density estimator is used. For "kde edge", a modified version of
#'     "kde" is used that allows for a possible point mass at the left edge of
#'     the marker distribution. For "histogram", a histogram estimator is used.
#' @param hist_bins If density_type="histogram", this value controls the number
#'     of bins used to construct the histogram.
#' @param dat The data object originally passed into \code{\link{est_ce}}. It is
#'     only necessary to pass this in if \code{density_type} is not set to
#'     "none".
#' @return A plot of CR/CVE estimates
#' @examples
#' data(hvtn505)
#' dat <- load_data(time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt",
#'                  marker="IgG_V2", covariates=c("age","BMI","bhvrisk"),
#'                  weights="wt", ph2="casecontrol", data=hvtn505)
#' \donttest{
#' ests_cox <- est_ce(dat=dat, type="Cox", t_0=578)
#' ests_np <- est_ce(dat=dat, type="NP", t_0=578)
#' plot_ce(ests_cox, ests_np)
#' }
#' @export
plot_ce <- function(..., which="CR", labels=NA, density_type="none",
                    hist_bins=NA, dat=NA) {

  if (!(which %in% c("CR", "CVE"))) {
    stop("`which` must equal one of c('CR','CVE').")
  }

  if (length(list(...))==0) {
    stop(paste0("One or more objects of class 'vaccine_est' must be passed int",
                "o `plot_ce`."))
  }

  if (density_type!="none" && !is.na(dat)) {
    stop("If `density_type` is not set to 'none', `dat` must also be provided.")
  }

  # !!!!! Add to examples:
  # plot_ce(ests_cox, ests_np, density=list(s=dat_v$s, weights=dat_v$weights))

  # To prevent R CMD CHECK notes
  x <- y <- ci_lower <- ci_upper <- curve <- NULL
  rm(x, y, ci_lower, ci_upper, curve)

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
    if (class(obj)!="vaccine_est") {
      stop(paste0("One or more of the objects passed into `plot_ce` is not of ",
                  "of class 'vaccine_est'."))
    }
    if (which=="CR") {
      if (class(obj$cr)=="NULL") {
        stop(paste0("CR estimates not present in one or more `vaccine_est` obj",
                    "ects; try rerunning `est_ce` with cr=TRUE."))
      } else {
        df_add <- data.frame(
          x = obj$cr$s,
          y = obj$cr$est,
          ci_lower = obj$cr$ci_lower,
          ci_upper = obj$cr$ci_upper,
          curve = labels[counter]
        )
      }
    } else if (which=="CVE") {
      if (class(obj$cve)=="NULL") {
        stop(paste0("CVE estimates not present in one or more `vaccine_est` ob",
                    "jects; try rerunning `est_ce` with cve=TRUE."))
      } else {
        df_add <- data.frame(
          x = obj$cve$s,
          y = obj$cve$est,
          ci_lower = obj$cve$ci_lower,
          ci_upper = obj$cve$ci_upper,
          curve = labels[counter]
        )
      }
    }
    df_plot <- rbind(df_plot,df_add)
    counter <- counter + 1
  }

  plot <- ggplot2::ggplot(df_plot, ggplot2::aes(x=x, y=y, color=curve, fill=curve))

  # if (!missing(density)) {
    # plot <- plot +
    #   ggplot2::geom_density(
    #     ggplot2::aes(x=x),
    #     data = data.frame(x=density$s),
    #     inherit.aes = F,
    #     fill = "forestgreen",
    #     alpha = 0.3,
    #     color = NA
    #   )
    # # plot
    # # inds <- which(!is.na(density$s))
    # # dens <- density(x=density$s[inds], weights=density$weights[inds])
    # # df_plot2 <- data.frame(x=dens$x, y=dens$y, color=NA, fill=NA)
    # # plot <- plot + geom_area(data=df_plot2, fill="forestgreen", alpha=0.5, color="white")
  # }

  # curve_colors <- c("darkgrey", "darkorchid3", "firebrick3", "deepskyblue3",
  #                   "darkgreen", "darkorange")
  curve_colors <- c("darkgreen", "darkorange")

  if (which=="CR") {
    labs <- list(title="Controlled Risk", y="Risk")
  } else {
    labs <- list(title="Controlled Vaccine Efficacy", y="CVE")
  }
  plot <- plot + ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=ci_lower, ymax=ci_upper), # color="darkorchid3", fill="darkorchid3"
                         alpha = 0.1, linetype = "dotted") +
    ggplot2::labs(title=labs$title, x="S", y=labs$y, color="",
                  fill="") +
    ggplot2::theme(
      legend.position = "bottom"
    )
  #

  # !!!!! Cut off at quantiles

  return(plot)

}



#' Trim data for plotting/reporting
#'
#' @description Removes a subset of estimates returned by \code{\link{est_ce}}
#' @param ests An object of class \code{"vaccine_est"} returned by
#'     \code{\link{est_ce}}.
#' @param dat The data object originally passed into \code{\link{est_ce}}.
#' @param quantiles A vector of length 2 representing the quantiles of the
#'     marker distribution at which to trim the data; if, for example,
#'     \code{quantiles=c(0.1,0.9)} is specified, values outside the 10% and 90%
#'     (weighted) quantiles of the marker distribution will be trimmed.
#' @return A modified copy of \code{ests} with the data trimmed.
#' @examples
#' data(hvtn505)
#' dat <- load_data(time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt",
#'                  marker="IgG_V2", covariates=c("age","BMI","bhvrisk"),
#'                  weights="wt", ph2="casecontrol", data=hvtn505)
#' \donttest{
#' ests_cox <- est_ce(dat=dat, type="Cox", t_0=578)
#' ests_cox <- trim(ests_cox, dat=dat, quantiles=c(0.1,0.9))
#' plot_ce(ests_cox, density_type="kde")
#' }
#' @export
trim <- function(ests, dat, quantiles) {

  dat_v <- dat[dat$a==1,]
  cutoffs <- quantile(dat_v$s, na.rm=T, probs=quantiles) # !!!!! Make this weighted

  for (i in c(1:length(ests))) {
    if (names(ests)[i] %in% c("cr", "cve")) {
      inds_to_keep <- which(ests[[i]]$s>=cutoffs[[1]] & ests[[i]]$s<=cutoffs[[2]])
      for (v in c("s", "est", "se", "ci_lower", "ci_upper")) {
        ests[[i]][[v]] <- ests[[i]][[v]][inds_to_keep]
      }
    }
  }

  return(ests)

}
