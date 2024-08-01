#' Plotting controlled effect curves
#'
#' @description Plot CR and/or CVE curves
#' @param ... One or more objects of class \code{"vaccine_est"} returned by
#'     \code{\link{est_ce}}.
#' @param which One of c("CR", "CVE"); controls whether to plot CR curves or CVE
#'     curves.
#' @param density_type One of c("none", "kde", "kde edge").
#'     Controls the type of estimator used for the background marker density
#'     plot. For "none", no density plot is displayed. For "kde", a weighted
#'     kernel density estimator is used. For "kde edge", a modified version of
#'     "kde" is used that allows for a possible point mass at the left edge of
#'     the marker distribution.
#' @param dat The data object originally passed into \code{\link{est_ce}}, used
#'     for plotting densities. It is only necessary to pass this in if
#'     \code{density_type} is not set to "none".
#' @param dat_alt Alternative data object; a list containing one or more
#'     dataframes, each of the form \code{data.frame(s=..., weights=...)}.
#'     Column \code{s} contains biomarker values and column \code{weights}
#'     contains corresponding two-phase sampling weights. This can be used as an
#'     alternative to specifying \code{dat}, and is particularly useful for
#'     plotting multiple densities on a single plot. If plotting multiple
#'     densities, the order of the dataframes should correspond to the order of
#'     \code{"vaccine_est"} objects passed in. See examples.
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
#' \donttest{
#' # Plot one curve
#' data(hvtn505)
#' dat <- load_data(time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt",
#'                  marker="IgG_V2", covariates=c("age","BMI","bhvrisk"),
#'                  weights="wt", ph2="casecontrol", data=hvtn505)
#' ests_cox <- est_ce(dat=dat, type="Cox", t_0=578)
#' plot_ce(ests_cox, density_type="kde", dat=dat)
#'
#' # Trim display of plot according to quantiles of the biomarker distribution
#' ests_cox_tr <- trim(ests_cox, dat=dat, quantiles=c(0.05,0.95))
#' plot_ce(ests_cox_tr, density_type="kde", dat=dat)
#'
#' # Plot multiple curves (same biomarker)
#' ests_np <- est_ce(dat=dat, type="NP", t_0=578)
#' plot_ce(ests_cox, ests_np, density_type="kde", dat=dat)
#'
#' # Plot multiple curves (two different biomarkers)
#' dat2 <- load_data(time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt",
#'                   marker="IgG_env", covariates=c("age","BMI","bhvrisk"),
#'                   weights="wt", ph2="casecontrol", data=hvtn505)
#' ests_cox2 <- est_ce(dat=dat2, type="Cox", t_0=578)
#' dat_alt <- list(
#'   data.frame(s=dat$s[dat$a==1], weights=dat$weights[dat$a==1]),
#'   data.frame(s=dat2$s[dat2$a==1], weights=dat2$weights[dat2$a==1])
#' )
#' plot_ce(ests_cox, ests_cox2, density_type="kde", dat_alt=dat_alt)
#' }
#' @export
plot_ce <- function(..., which="CR", density_type="none", dat=NA, dat_alt=NA,
                    zoom_x="zoom in", zoom_y="zoom out") {

  # Error handling
  if (!(which %in% c("CR", "CVE"))) {
    stop("`which` must equal one of c('CR','CVE').")
  }
  if (length(list(...))==0) {
    stop(paste0("One or more objects of class 'vaccine_est' must be passed int",
                "o `plot_ce`."))
  }
  if (density_type!="none" && missing(dat) && missing(dat_alt)) {
    stop(paste0("If `density_type` is not set to 'none', either `dat` or `dat_",
                "alt` must also be provided."))
  }
  if (!missing(dat) && !missing(dat_alt)) {
    stop("Specify either `dat` or `dat_alt`, but not both.")
  }
  if (!missing(dat)) {
    if (!methods::is(dat,"vaccine_dat")) {
      stop(paste0("`dat` must be an object of class 'vaccine_dat' returned by ",
                  "load_data()."))
    }
  }
  if (!missing(dat_alt)) {
    violations <- !methods::is(dat_alt,"list") ||
      !methods::is(dat_alt[[1]],"data.frame") ||
      !identical(sort(names(dat_alt[[1]])), c("s", "weights"))
    if (violations) {
      stop(paste0("`dat_alt` must be a list containing one or more dataframes,",
                  " of the form `list(data.frame(s=..., weights=...))`."))
    }
    use_dat_alt <- TRUE
  } else {
    use_dat_alt <- FALSE
  }
  if (!(length(zoom_x) %in% c(1,2)) ||
      (length(zoom_x)==2 && !methods::is(zoom_x, "numeric")) ||
      (length(zoom_x)==1 && !(zoom_x %in% c("zoom out", "zoom in")))) {
    stop(paste0("`zoom_x` must equal either 'zoom in' or 'zoom out', or be a n",
                "umeric vector of length 2."))
  }
  if (length(zoom_x)==1 && zoom_x=="zoom out" && missing(dat)) {
    stop("If `zoom_x` is set to 'zoom out', `dat` must be provided as well.")
  }
  if (!(length(zoom_y) %in% c(1,2)) ||
      (length(zoom_y)==2 && !methods::is(zoom_y, "numeric")) ||
      (length(zoom_y)==1 && zoom_y!="zoom out")) {
    stop(paste0("`zoom_y` must equal 'zoom out' or be a numeric vector ",
                "of length 2."))
  }

  # To prevent R CMD CHECK notes
  x <- y <- ci_lower <- ci_upper <- curve <- ymin <- ymax <- NULL
  rm(x, y, ci_lower, ci_upper, curve, ymin, ymax)

  df_plot <- as_table(..., which=which)

  # Color scheme
  curve_colors <- c("deepskyblue3", "darkorchid3", "darkgreen", "darkorange",
                    "firebrick3", "darkgrey")

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

  # Construct KDE dataframe
  if (density_type!="none") {

    if (!is.null(attr(list(...)[[1]], "cr_placebo_arm")) && !use_dat_alt) {
      dat$a <- 1-dat$a
    }

    if (!use_dat_alt) { dat_alt <- list(dat[dat$a==1,]) }

    kde_area <- c()

    for (i in c(1:length(dat_alt))) {

      dat_ <- dat_alt[[i]]
      min_s <- min(dat_$s, na.rm=T)
      p_edge <- sum(dat_$weights*In(dat_$s==min_s), na.rm=T)/nrow(dat_)
      if (p_edge<0.03 & density_type=="kde edge") { density_type <- "kde" }
      dens_height <- 0.6 * (z_y[2]-z_y[1])

      if (density_type=="kde") {

        df_dens <- data.frame(
          s = dat_$s[!is.na(dat_$s)],
          weights = dat_$weights[!is.na(dat_$s)]
        )
        df_dens$weights <- df_dens$weights / sum(df_dens$weights)
        dens <- suppressWarnings(stats::density(
          x = df_dens$s,
          bw = "ucv",
          weights = df_dens$weights
        ))

      } else {

        df_dens <- data.frame(
          s = dat_$s[!is.na(dat_$s) & dat_$s!=min_s],
          weights = dat_$weights[!is.na(dat_$s) & dat_$s!=min_s]
        )
        df_dens$weights <- df_dens$weights / sum(df_dens$weights)
        dens <- suppressWarnings(stats::density(
          x = df_dens$s,
          bw = "ucv",
          weights = df_dens$weights
        ))
        dens$y <- dens$y * (1-p_edge)

        plot_width <- z_x[2]-z_x[1]
        rect_x <- c(min_s-0.025*plot_width, min_s+0.025*plot_width)
        rect_y <- p_edge / (rect_x[2]-rect_x[1])
        inds_to_remove <- dens$x>rect_x[2]
        dens$x <- dens$x[inds_to_remove]
        dens$y <- dens$y[inds_to_remove]
        dens$x[length(dens$x)+1] <- rect_x[1]
        dens$y[length(dens$y)+1] <- rect_y
        dens$x[length(dens$x)+1] <- rect_x[2]
        dens$y[length(dens$y)+1] <- rect_y
        dens$x[length(dens$x)+1] <- rect_x[2] + plot_width/10^5
        dens$y[length(dens$y)+1] <- z_y[1]

      }

      kde_data <- data.frame(
        x = dens$x,
        ymin = z_y[1],
        ymax = dens_height * (dens$y/max(dens$y)) + z_y[1]
      )
      assign(x=paste0("kde_data_",i), value=kde_data)

      inds <- c(2:nrow(kde_data))
      kde_area_x <- kde_data$x[inds] - kde_data$x[inds-1]
      kde_area_y <- kde_data$ymax[inds] - kde_data$ymin[inds]
      kde_area[i] <- sum(kde_area_x*kde_area_y)

    }

    for (i in c(1:length(dat_alt))) {

      # If there are multiple densities, scale them so that they have the same area
      kde_data <- get(paste0("kde_data_",i))
      area_mult <- min(kde_area)/kde_area[i]
      kde_data$ymax <- area_mult*(kde_data$ymax-kde_data$ymin) + kde_data$ymin

      # Plot background KDE
      fill_color <- ifelse(length(dat_alt)==1, "orange", curve_colors[i])
      plot <- plot + ggplot2::geom_ribbon(
        ggplot2::aes(x=x, ymin=ymin, ymax=ymax),
        data = kde_data,
        inherit.aes = F,
        color = "white",
        fill = fill_color,
        alpha = 0.2
      )

    }

  }

  # Primary plot
  plot <- plot +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin=ci_lower, ymax=ci_upper), # color="darkorchid3", fill="darkorchid3"
                         alpha = 0.05, linetype = "dotted") +
    ggplot2::labs(title=labs$title, x="S", y=labs$y, color="",
                  fill="") +
    ggplot2::scale_color_manual(values=curve_colors) +
    ggplot2::scale_fill_manual(values=curve_colors) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color="#bbbbbb", fill=NA)
    )

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
#' \donttest{
#' data(hvtn505)
#' dat <- load_data(time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt",
#'                  marker="IgG_V2", covariates=c("age","BMI","bhvrisk"),
#'                  weights="wt", ph2="casecontrol", data=hvtn505)
#' ests_cox <- est_ce(dat=dat, type="Cox", t_0=578)
#' ests_cox_tr <- trim(ests_cox, dat=dat, quantiles=c(0.05,0.95))
#' plot_ce(ests_cox_tr, density_type="kde", dat=dat)
#' }
#' @export
trim <- function(ests, dat, quantiles) {

  dat_v <- dat[dat$a==1,]
  cutoffs <- stats::quantile(dat_v$s, na.rm=T, probs=quantiles) # !!!!! Make this weighted quantile

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



#' Create table of estimates
#'
#' @description Format estimates returned by \code{\link{est_ce}} as a table
#' @param ... One or more objects of class \code{"vaccine_est"} returned by
#'     \code{\link{est_ce}}.
#' @param which One of c("CR", "CVE"); controls whether the table contains CR or
#'     CVE values.
#' @return A table of CR or CVE values
#' @examples
#' data(hvtn505)
#' dat <- load_data(time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt",
#'                  marker="IgG_V2", covariates=c("age","BMI","bhvrisk"),
#'                  weights="wt", ph2="casecontrol", data=hvtn505)
#' \donttest{
#' ests_cox <- est_ce(dat=dat, type="Cox", t_0=578)
#' ests_np <- est_ce(dat=dat, type="NP", t_0=578)
#' ests_table <- as_table(ests_cox, ests_np)
#' head(ests_table)
#' }
#' @export
as_table <- function(..., which="CR") {

  # !!!!! Move all error handling into a function

  df_ests <- data.frame(
    x = double(),
    y = double(),
    ci_lower = double(),
    ci_upper = double(),
    curve = character()
  )

  labels <- sapply(substitute(list(...))[-1], deparse)

  counter <- 1
  for (obj in list(...)) {
    if (!methods::is(obj, "vaccine_est")) {
      stop(paste0("One or more of the objects passed into `plot_ce` is not of ",
                  "of class 'vaccine_est'."))
    }
    if (which=="CR") {
      if (methods::is(obj$cr, "NULL")) {
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
      if (methods::is(obj$cve, "NULL")) {
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
    df_ests <- rbind(df_ests,df_add)
    counter <- counter + 1
  }

  return(df_ests)

}
