#' Plotting controlled effect curves
#'
#' @description Plot CR and/or CVE curves
#' @param ... One or more objects of class \code{"vaccine_est"} returned by
#'     \code{\link{est_ce}}.
#' @param which One of c("CR", "CVE"); controls whether to plot CR curves or CVE
#'     curves.
#' @param labels A character vector of labels corresponding to the estimate
#'     objects.
#' @param density_type One of c("none", "kde", "kde edge").
#'     Controls the type of estimator used for the background marker density
#'     plot. For "none", no density plot is displayed. For "kde", a weighted
#'     kernel density estimator is used. For "kde edge", a modified version of
#'     "kde" is used that allows for a possible point mass at the left edge of
#'     the marker distribution.
#' @param dat The data object originally passed into \code{\link{est_ce}}. It is
#'     only necessary to pass this in if \code{density_type} is not set to
#'     "none".
#' @param zoom_x Either one of c("zoom in", "zoom out") or a vector of
#'     length 2. Controls the zooming on the X-axis. The default "zoom in" will
#'     set the zoom limits to the plot estimates. Choosing "zoom out" will set
#'     the zoom limits to show the entire distribution of the marker. Entering a
#'     vector of length 2 will set the left and right zoom limits explicitly.
#' @param zoom_y Either "zoom out" or a vector of length 2. Controls the zooming
#'     on the Y-axis. The default "zoom out" will show the entire vertical range
#'     of the estimates. Entering a vector of length 2 will set the lower and
#'     upper zoom limits explicitly.
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
plot_ce <- function(..., which="CR", labels=NA, density_type="none", dat=NA,
                    zoom_x="zoom in", zoom_y="zoom out") {

  # !!!!! Check dat
  # !!!!! Move all error handling into a function

  # Error handling
  if (!(which %in% c("CR", "CVE"))) {
    stop("`which` must equal one of c('CR','CVE').")
  }
  if (length(list(...))==0) {
    stop(paste0("One or more objects of class 'vaccine_est' must be passed int",
                "o `plot_ce`."))
  }
  if (density_type!="none" && missing(dat)) {
    stop("If `density_type` is not set to 'none', `dat` must also be provided.")
  }
  if (!(length(zoom_x) %in% c(1,2)) ||
      (length(zoom_x)==2 && !(class(zoom_x)=="numeric")) ||
      (length(zoom_x)==1 && !(zoom_x %in% c("zoom out", "zoom in")))) {
    stop(paste0("`zoom_x` must equal either 'zoom in' or 'zoom out', or be a n",
                "umeric vector of length 2."))
  }
  if (length(zoom_x)==1 && zoom_x=="zoom out" && missing(dat)) {
    stop("If `zoom_x` is set to 'zoom out', `dat` must be provided as well.")
  }
  if (!(length(zoom_y) %in% c(1,2)) ||
      (length(zoom_y)==2 && !(class(zoom_y)=="numeric")) ||
      (length(zoom_y)==1 && zoom_y!="zoom out")) {
    stop(paste0("`zoom_y` must equal 'zoom out' or be a numeric vector ",
                "of length 2."))
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
    ci_upper = double(),
    curve = character()
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

  # Set X-axis zoom values
  z_x <- rep(NA,2)
  if (length(zoom_x)==2) {
    z_x <- zoom_x
  } else {
    if (zoom_x=="zoom in") {
      z_x[1] <- min(df_plot$x)
      z_x[2] <- max(df_plot$x)
    } else if (zoom_x=="zoom out") {
      dat_v <- dat[dat$a==1,]
      z_x[1] <- min(dat_v$s, na.rm=T)
      z_x[2] <- max(dat_v$s, na.rm=T)
    }
    # Add 5% padding to zoom
    z_x[1] <- z_x[1] - 0.05*(z_x[2]-z_x[1])
    z_x[2] <- z_x[2] + 0.05*(z_x[2]-z_x[1])
  }

  # Set Y-axis zoom values
  z_y <- rep(NA,2)
  if (length(zoom_y)==2) {
    z_y <- zoom_y
  } else if (zoom_y=="zoom out") {
    zz <- dplyr::filter(df_plot, x>=z_x[1] & x<=z_x[2])
    z_y[1] <- ifelse(which=="CVE", min(zz$ci_lower, na.rm=T), 0)
    z_y[2] <- max(zz$ci_upper, na.rm=T)
    # Add 5% padding to zoom
    z_y[1] <- z_y[1] - 0.05*(z_y[2]-z_y[1])
    z_y[2] <- z_y[2] + 0.05*(z_y[2]-z_y[1])
  }

  # Construct KDE dataframe
  if (density_type!="none") {

    dat_v <- dat[dat$a==1,]
    min_s <- min(dat_v$s, na.rm=T)
    p_edge <- mean(dat_v$s==min_s, na.rm=T) # !!!!! Make this weighted
    if (p_edge<0.03 & density_type=="kde edge") { density_type <- "kde" }

    if (density_type=="kde") {

      dens_height <- 0.6 * (z_y[2]-z_y[1])
      df_dens <- data.frame(
        s = dat_v$s[!is.na(dat_v$s)],
        weights = dat_v$weights[!is.na(dat_v$s)]
      )
      df_dens$weights <- df_dens$weights / sum(df_dens$weights)
      dens <- suppressWarnings(stats::density(
        x = df_dens$s,
        bw = "ucv",
        # adjust = 2, # !!!!!
        weights = df_dens$weights
      ))
      kde_data <- data.frame(
        x = dens$x,
        ymin = z_y[1],
        ymax = dens_height * (dens$y/max(dens$y)) + z_y[1]
      )

    } else {

      stop("TO DO")

    }

  }

  # Set plot labels
  if (which=="CR") {
    labs <- list(title="Controlled Risk", y="Risk")
  } else {
    labs <- list(title="Controlled Vaccine Efficacy", y="CVE")
  }

  # Set up ggplot2 object
  plot <- ggplot2::ggplot(
    df_plot,
    ggplot2::aes(x=x, y=y, color=curve, fill=curve)
  ) +
    ggplot2::coord_cartesian(xlim=z_x, ylim=z_y, expand=F)

  # Plot background KDE
  if (density_type!="none") {
    plot <- plot + ggplot2::geom_ribbon(
      ggplot2::aes(x=x, ymin=ymin, ymax=ymax),
      data = kde_data,
      inherit.aes = F,
      color = "white",
      fill = "orange", # "forestgreen"
      alpha = 0.3
    )
  }

  # Primary plot
  plot <- plot +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=ci_lower, ymax=ci_upper), # color="darkorchid3", fill="darkorchid3"
                         alpha = 0.1, linetype = "dotted") +
    ggplot2::labs(title=labs$title, x="S", y=labs$y, color="",
                  fill="") +
    ggplot2::theme(
      legend.position = "bottom",
      panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color="#bbbbbb", fill=NA)
    )

  # Implement these eventually
  curve_colors <- c("darkgrey", "darkorchid3", "firebrick3", "deepskyblue3",
                    "darkgreen", "darkorange")
  # curve_colors <- c("darkgreen", "darkorange")

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
