#' Load and format data object
#'
#' @description This function takes in user-supplied data and returns a data
#'     object that can be read in by \code{est_np()}, \code{est_cox()},
#'     and \code{hyptest_np()}. Data is expected to come from a case-cohort
#'     sampling design and be filtered to include all phase-one individuals.
#' @param time A numeric vector of observed event or censoring times.
#' @param event A vector of binary values corresponding to whether the observed
#'     time represents an event time (1) or a censoring time (0). Accepts either
#'     integer (0/1) or Boolean (T/F) values.
#' @param vacc A vector of binary values corresponding to whether the individual
#'     is in the vaccine group (1) or the placebo group (0). Accepts either
#'     integer (0/1) or Boolean (T/F) values.
#' @param marker A numeric vector of biomarker values.
#' @param covariates A dataframe. Columns values should be either numeric,
#'     binary, or factors. Character columns will be converted into factors.
#' @param weights A numeric vector of inverse-probability-of-sampling (IPS)
#'     weights.
#' @param ph2 A vector of binary values corresponding to whether the individual
#'     is in the phase-two cohort (1) or not (0). Accepts either integer (0/1)
#'     or Boolean (T/F) values.
#' @return A list containing the following: \itemize{
#'     \item{\code{one}: asdf}
#'     \item{\code{two}: asdf}
#'     \item{\code{three}: asdf}
#' }
#' @examples
#' print("to do")
#' @export
load_data <- function(
  time, event, vacc, marker, covariates, weights, ph2
) {

  # !!!!! Add strata; this should be converted to an integer such that unique strata are numbered 1:(# strata)

  # Input checks: time
  n <- length(time)
  if (n==0) { stop("`time` vector cannot be of length zero.") }
  if (any(is.na(event))) {
    stop("NA values not allowed in `time` vector.")
  }
  if (!is.numeric(event)) { stop("`time` vector must be numeric.") }

  # Input checks: event
  if (length(event)!=n) {
    stop("`time` and `event` vectors must be of the same length.")
  }
  if (any(is.na(event))) {
    stop("NA values not allowed in `event` vector.")
  }
  if (any(!(event %in% c(0,1,F,T)))) {
    stop("`time` vector must only contain binary values (either F/T or 0/1).")
  }
  event <- as.integer(event)

  # Input checks: vacc
  if (length(vacc)!=n) {
    stop("`time` and `vacc` vectors must be of the same length.")
  }
  if (any(is.na(vacc))) {
    stop("NA values not allowed in `vacc` vector.")
  }
  if (any(!(vacc %in% c(0,1,F,T)))) {
    stop("`time` vector must only contain binary values (either F/T or 0/1).")
  }
  vacc <- as.integer(vacc)
  .groups <- ifelse(any(vacc==0) && any(vacc==1), "both",
                    ifelse(any(vacc==1), "vaccine", "placebo"))

  # Input checks: marker
  if (length(marker)!=n) {
    stop("`time` and `event` vectors must be of the same length.")
  }
  if (!is.numeric(marker)) { stop("`marker` vector must be numeric.") }

  # Input checks: covariates
  if (class(covariates)!="data.frame") {
    stop("`covariates` must be a data frame.")
  }
  if (nrow(covariates)!=n) {
    stop(paste0("Length of `time` vector must equal the number of rows in the ",
                "`covariates` dataframe.")) }
  if (any(is.na(covariates))) {
    stop("NA values not allowed in `covariates` dataframe.")
  }
  names(covariates) <- paste0("x", c(1:length(covariates)))
  # !!!!! convert factors to dummy columns

  # Input checks: weights
  if (length(weights)!=n) {
    stop("`time` and `weights` vectors must be of the same length.")
  }
  if (!is.numeric(weights)) { stop("`weights` vector must be numeric.") }

  # Input checks: ph2
  if (length(ph2)!=n) {
    stop("`time` and `ph2` vectors must be of the same length.")
  }
  if (any(is.na(ph2))) {
    stop("NA values not allowed in `ph2` vector")
  }
  if (any(!(ph2 %in% c(0,1,F,T)))) {
    stop("`ph2` vector must only contain binary values (either F/T or 0/1).")
  }
  ph2 <- as.integer(ph2)

  if (.groups %in% c("vaccine", "both")) {

    # Create data object
    .ind_v <- which(vacc==1)
    df_vc <- list(
      "y" = time[.ind_v],
      "delta" = event[.ind_v],
      "s" = marker[.ind_v],
      "x" = covariates[.ind_v,],
      "weights" = ph2*weights[.ind_v],
      "strata" = factor(ph2*weights[.ind_v]),
      "z" = ph2[.ind_v]
    )
    attr(df_vc, "n_orig") <- length(df_vc$z)
    attr(df_vc, "dim_x") <- length(covariates)

    # Stabilize weights (rescale to sum to sample size)
    .stb_v <- sum(df_vc$weights) / length(df_vc$z)
    df_vc$weights <- df_vc$weights / .stb_v

  } else {
    df_vc <- NA
  }

  if (.groups %in% c("placebo", "both")) {

    # Create data object
    .ind_p <- which(vacc==1)
    df_pl <- list(
      "y" = time[.ind_p],
      "delta" = event[.ind_p],
      "s" = marker[.ind_p],
      "x" = covariates[.ind_p,],
      "weights" = ph2*weights[.ind_p],
      "strata" = factor(ph2*weights[.ind_p]),
      "z" = ph2[.ind_p]
    )
    attr(df_vc, "n_orig") <- length(df_pl$z)

    # Stabilize weights (rescale to sum to sample size)
    .stb_p <- sum(df_pl$weights) / length(df_pl$z)
    df_pl$weights <- df_pl$weights / .stb_p

  } else {
    df_pl <- NA
  }

  # Create and return data object
  dat <- list("v"=df_vc, "p"=df_pl)
  class(dat) <- "dat_vaccine"
  return(dat)

}
