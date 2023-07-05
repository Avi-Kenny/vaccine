##' Load and format data object
#'
#' @description This function takes in user-supplied data and returns a data
#'     object that can be read in by \code{\link{summary_stats}},
#'     \code{\link{est_ce}}, and other estimation functions. Data is expected to
#'     come from a two-phase sampling design and should be filtered to include
#'     all phase-one individuals.
#' @param time A character string; the name of the numeric variable representing
#'     observed event or censoring times.
#' @param event A character string; the name of the binary variable
#'     corresponding to whether the observed time represents an event time (1)
#'     or a censoring time (0). Either integer (0/1) or Boolean (T/F) values are
#'     allowed.
#' @param vacc A character string; the name of the binary variable denoting
#'     whether the individual is in the vaccine group (1) or the placebo group
#'     (0). Accepts either integer (0/1) or Boolean (T/F) values.
#' @param marker A character string; the name of the numeric variable of
#'     biomarker values.
#' @param covariates A character vector; the names of the covariate columns.
#'     Columns values should be either numeric, binary, or factors. Character
#'     columns will be converted into factors.
#' @param weights A character string; the name of the numeric variable
#'     containing inverse-probability-of-sampling (IPS) weights.
#' @param ph2 A character string; the name of the binary variable representing
#'     whether the individual is in the phase-two cohort (1) or not (0). Accepts
#'     either integer (0/1) or Boolean (T/F) values.
#' @param strata A character string; the name of the variable containing strata
#'     identifiers (for two-phase sampling strata).
#' @param data A dataframe containing the vaccine trial data.
#' @return An object of class \code{vaccine_dat}.
#' @examples
#' \dontrun{
#' data(hvtn505)
#' dat <- load_data(time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt",
#'                  marker="logpctpos_scaled", covariates=c("age","BMI","bhvrisk"),
#'                  weights="wt", ph2="casecontrol", data=hvtn505)
#' }
#' @export
load_data <- function(
  time, event, vacc, marker, covariates, weights, ph2, strata=NA, data
) {

  # To prevent R CMD CHECK notes
  .time <- .marker <- .covariates <- NULL; rm(.time,.marker,.covariates);

  # Input validation
  {

    if (!methods::is(data,"data.frame")) { stop("`data` must be a data frame.") }
    if (nrow(data)==0) { stop("`data` is an empty data frame.") }

    for (arg in c("time", "event", "vacc", "marker", "covariates", "weights",
                  "ph2", "strata")) {

      var <- get(arg)
      if (!(arg=="strata" && missing(strata))) {

        # Variable is a character string specifying variable(s) in `data`
        is_string <- methods::is(var,"character")
        length_one <- as.logical(length(var)==1)
        in_df <- all(as.logical(var %in% names(data)))
        if (arg=="covariates") {
          if (!(is_string&&in_df)) {
            stop(paste0("`", arg, "` must be a vector of character strings spe",
                        "cifying one or more variables in `data`."))
          }
        } else {
          if (!(is_string&&length_one&&in_df)) {
            stop(paste0("`", arg, "` must be a character string specifying a s",
                        "ingle variable in `data`."))
          }
        }

        # Assign column(s) to val
        if (arg=="covariates") {
          val <- data[,var, drop=F]
        } else {
          val <- data[,var]
        }

        # No missing values allowed (except marker)
        if (!(arg %in% c("marker", "weights"))) {
          if (any(is.na(val))) { stop("NA values not allowed in `", arg, "`.") }
        }

        # Validate: `time`, `marker`, `weights`
        if (arg %in% c("time", "marker", "weights")) {
          if (!is.numeric(val)) { stop(paste0("`", arg, "` must be numeric.")) }
        }

        # Validate: `event`, `vacc`, `ph2`
        if (arg %in% c("event", "vacc", "ph2")) {
          if (any(!(val %in% c(0,1,F,T)))) {
            stop(paste0("`", arg, "` must only contain binary values (either T",
                        "/F or 1/0)."))
          }
        }

      }

      assign(x=paste0(".",arg), value=val)

    }

  }

  # Convert binary variables to integers (if specified as boolean)
  .event <- as.integer(.event)
  .vacc <- as.integer(.vacc)
  .ph2 <- as.integer(.ph2)

  .groups <- ifelse(any(.vacc==0) && any(.vacc==1), "both",
                    ifelse(any(.vacc==1), "vaccine", "placebo"))

  .weights <- ifelse(is.na(.weights), 0, .weights)

  # !!!!! convert factors to dummy columns; import code from VaxCurve
  # !!!!! Check what processing we have to do to strata (e.g. convert to integers)
  # !!!!! Store two copies of covariates; one for Cox model and one for NPCVE etc.
  # !!!!! Also maybe store a combined version of the dataset (or have a helper function to combine)?
  # Change n_orig to n or n_ph1

  if (.groups %in% c("vaccine", "both")) {

    # Create strata (if not given)
    .ind_v <- which(.vacc==1)
    if(is.na(.strata[[1]])) {
      .strata <- as.integer(factor(.weights[.ind_v]))
    } else {
      .strata <- as.integer(factor(.strata[.ind_v]))
    }

    # Create data object
    df_vc <- list(
      "y" = .time[.ind_v],
      "delta" = .event[.ind_v],
      "s" = .marker[.ind_v],
      "x" = .covariates[.ind_v,, drop=F],
      "weights" = .ph2[.ind_v]*.weights[.ind_v],
      "strata" = .strata,
      "z" = .ph2[.ind_v]
    )
    names(df_vc$x) <- paste0("x", c(1:length(df_vc$x)))
    attr(df_vc, "n_orig") <- length(df_vc$z)
    attr(df_vc, "dim_x") <- length(.covariates)

    # Stabilize weights (rescale to sum to sample size)
    .stb_v <- sum(df_vc$weights) / length(df_vc$z)
    df_vc$weights <- df_vc$weights / .stb_v

  } else {
    df_vc <- list()
  }

  if (.groups %in% c("placebo", "both")) {

    # Create strata (if not given)
    .ind_p <- which(.vacc==0)
    if(is.na(.strata[[1]])) {
      .strata <- as.integer(factor(.weights[.ind_p]))
    } else {
      .strata <- as.integer(factor(.strata[.ind_p]))
    }

    # Create data object
    df_pl <- list(
      "y" = .time[.ind_p],
      "delta" = .event[.ind_p],
      "s" = .marker[.ind_p],
      "x" = .covariates[.ind_p,, drop=F],
      "weights" = .ph2[.ind_p]*.weights[.ind_p],
      "strata" = .strata,
      "z" = .ph2[.ind_p]
    )
    names(df_pl$x) <- paste0("x", c(1:length(df_pl$x)))
    attr(df_pl, "n_orig") <- length(df_pl$z)

    # Stabilize weights (rescale to sum to sample size)
    .stb_p <- sum(df_pl$weights) / length(df_pl$z)
    df_pl$weights <- df_pl$weights / .stb_p

  } else {
    df_pl <- list()
  }

  # Create and return data object
  dat <- list("v"=df_vc, "p"=df_pl)
  class(dat) <- "vaccine_dat"
  attr(dat, "groups") <- .groups
  attr(dat, "covariate_names") <- covariates
  return(dat)

}
