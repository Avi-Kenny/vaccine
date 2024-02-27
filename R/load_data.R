##' Load and format data object
#'
#' @description This function takes in user-supplied data and returns a data
#'     object that can be read in by \code{\link{summary_stats}},
#'     \code{\link{est_ce}}, \code{\link{est_med}}, and other estimation
#'     functions. Data is expected to come from a vaccine clinical trial,
#'     possibly involving two-phase sampling and possibly including a biomarker
#'     of interest.
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
#' @param covariates_ph2 A boolean; if at least one of the covariates is
#'     measured only in the phase-two cohort, set this to TRUE.
#' @return An object of class \code{vaccine_dat}.
#' @examples
#' data(hvtn505)
#' dat <- load_data(time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt",
#'                  marker="IgG_V2", covariates=c("age","BMI","bhvrisk"),
#'                  weights="wt", ph2="casecontrol", data=hvtn505)
#' @export
load_data <- function(
  time, event, vacc, marker, covariates, weights, ph2, strata=NA, data,
  covariates_ph2=FALSE
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
      # if (!(arg=="strata")) { # For testing

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
          val2 <- list()

          # Convert factor/character columns
          for (col in names(val)) {
            if (is.numeric(val[,col])) {
              val2[[col]] <- val[,col]
            } else if (is.logical(val[,col])) {
              val2[[col]] <- as.integer(val[,col])
            } else if (is.factor(val[,col]) || is.character(val[,col])) {
              tmp_col <- as.integer(as.factor(val[,col]))
              tmp_unique <- unique(tmp_col)
              tmp_size <- length(tmp_unique)
              for (i in c(1:(tmp_size-1))) {
                col_new <- paste0(col, "_", i)
                val2[[col_new]] <- as.integer(tmp_col==tmp_unique[i])
              }
            } else {
              stop(paste0("The data type of column `", col, "` is not supported."))
            }
          }
          val <- as.data.frame(val2)
          x_names <- names(val)
        } else {
          val <- data[,var]
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

        # No missing values allowed (except marker)
        if (!(arg %in% c("marker", "weights", "covariates"))) {
          if (any(is.na(val))) { stop("NA values not allowed in `", arg, "`.") }
        } else if (arg=="covariates") {
          cond_1 <- any(is.na(val))
          cond_2 <- any(is.na(val[as.logical(data[,ph2]),]))
          msg_1 <- "NA values not allowed in `covariates` (if covariates_ph2==F)."
          msg_2 <- paste0("NA values only allowed in `covariates` for which ph2==F (if covariates_ph2==T).")
          if (cond_1 && !covariates_ph2) { stop(msg_1) }
          if (cond_2 && covariates_ph2) { stop(msg_2) }
        }

        assign(x=paste0(".",arg), value=val)

      } else {

        .strata <- NA

      }

    }

  }

  # Convert binary variables to integers (if specified as boolean)
  .event <- as.integer(.event)
  .vacc <- as.integer(.vacc)
  .ph2 <- as.integer(.ph2)

  # Create additional variables
  .groups <- ifelse(any(.vacc==0) && any(.vacc==1), "both",
                    ifelse(any(.vacc==1), "vaccine", "placebo"))

  .weights <- ifelse(is.na(.weights), 0, .weights)
  .dim_x <- length(.covariates)
  .n_v <- sum(.vacc)
  .n_p <- sum(1-.vacc)

  # Rename covariate dataframe to c("x1", "x2", ...)
  names(.covariates) <- paste0("x", c(1:.dim_x))

  # !!!!! Warning/error if there is only one unique level or too many (>15) unique levels
  # !!!!! Warning/error if there are numbers passed in as character strings
  # !!!!! Check what processing we have to do to strata (e.g. convert to integers)
  # !!!!! Convert booleans to integers (0 or 1) for all coolumns (incl covariates)
  # !!!!! Store two copies of covariates; one for Cox model and one for NPCVE etc.
  # !!!!! Also maybe store a combined version of the dataset (or have a helper function to combine)?

  if (.groups %in% c("vaccine", "both")) {

    # Create strata (if not given)
    .ind_v <- which(.vacc==1)
    if(is.na(.strata[[1]])) {
      .strata <- as.integer(factor(.weights[.ind_v]))
    } else {
      .strata <- as.integer(factor(.strata[.ind_v]))
    }

    # Create data object
    df_v <- cbind(
      .covariates[.ind_v,, drop=F], "y"=.time[.ind_v], "delta"=.event[.ind_v],
      "s"=.marker[.ind_v], "weights"=.ph2[.ind_v]*.weights[.ind_v],
      "strata"=.strata, "z"=.ph2[.ind_v], "a"=1
    )

    # Stabilize weights (rescale to sum to sample size)
    .stb_v <- sum(df_v$weights) / .n_v
    df_v$weights <- df_v$weights / .stb_v

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
    df_p <- cbind(
      .covariates[.ind_p,, drop=F], "y"=.time[.ind_p], "delta"=.event[.ind_p],
      "s"=.marker[.ind_p], "weights"=.ph2[.ind_p]*.weights[.ind_p],
      "strata"=.strata, "z"=.ph2[.ind_p], "a"=0
    )

    # Stabilize weights (rescale to sum to sample size)
    # !!!!! Consider making weights NA in placebo group
    .stb_p <- sum(df_p$weights) / .n_p
    df_p$weights <- df_p$weights / .stb_p

  }

  if (.groups=="vaccine") {
    dat <- df_v
  } else if (.groups=="placebo") {
    dat <- df_p
  } else if (.groups=="both") {
    dat <- rbind(df_v, df_p)
  }

  # Create and return data object
  class(dat) <- c("data.frame", "vaccine_dat")
  attr(dat, "groups") <- .groups
  attr(dat, "covariate_names") <- x_names
  attr(dat, "covariates_ph2") <- covariates_ph2
  attr(dat, "dim_x") <- .dim_x
  attr(dat, "n") <- .n_v+.n_p
  attr(dat, "n_vacc") <- .n_v
  attr(dat, "n_vacc2") <- sum(df_v$z)
  attr(dat, "n_plac") <- .n_p
  return(dat)

}
