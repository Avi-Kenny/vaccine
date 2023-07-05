#' Calculate summary statistics
#'
#' @description TO DO
#' @param dat A data object returned by `load_data`.
#' @param quietly Boolean. If true, output will not be printed.
#' @return A list containing values of various summary statistics.
#' @examples
#' \dontrun{
#' data(hvtn505)
#' dat <- load_data(time="HIVwk28preunblfu", event="HIVwk28preunbl", vacc="trt",
#'                  marker="IgG_V2", covariates=c("age","BMI","bhvrisk"),
#'                  weights="wt", ph2="casecontrol", data=hvtn505)
#' summary_stats(dat=dat)
#' }
#' @export
summary_stats <- function(dat, quietly=F) {

  if (!methods::is(dat,"vaccine_dat")) {
    stop(paste0("`dat` must be an object of class 'vaccine_dat' returned by lo",
                "ad_data()."))
  }

  # !!!!! Checks based on identical(dat$v,list()) and identical(dat$p,list())

  num_ph1_subj_v <- length(dat$v$delta)
  num_ph1_subj_p <- length(dat$p$delta)
  num_ph2_subj_v <- sum(dat$v$z)
  num_ph2_subj_p <- sum(dat$p$z)
  num_ph1_events_v <- sum(dat$v$delta==1)
  num_ph1_events_p <- sum(dat$p$delta==1)
  num_ph2_events_v <- sum(dat$v$z==1 & dat$v$delta==1)
  num_ph2_events_p <- sum(dat$p$z==1 & dat$p$delta==1)

  ss <- list(
    "num_ph1_subj_v" = num_ph1_subj_v,
    "num_ph1_subj_p" = num_ph1_subj_p,
    "num_ph2_subj_v" = num_ph2_subj_v,
    "num_ph2_subj_p" = num_ph2_subj_p,
    "num_ph1_events_v" = num_ph1_events_v,
    "num_ph1_events_p" = num_ph1_events_p,
    "num_ph2_events_v" = num_ph2_events_v,
    "num_ph2_events_p" = num_ph2_events_p,
    "prop_ph1_events_v" = round(num_ph1_events_v/num_ph1_subj_v, 5),
    "prop_ph1_events_p" = round(num_ph1_events_p/num_ph1_subj_p, 5),
    "prop_ph2_events_v" = round(num_ph2_events_v/num_ph2_subj_v, 5),
    "prop_ph2_events_p" = round(num_ph2_events_p/num_ph2_subj_p, 5)
  )

  if (!quietly) {

    # Number of subjects
    message(paste("Number of subjects (vaccine group, phase-1):",
                  ss$num_ph1_subj_v))
    message(paste("Number of subjects (placebo group, phase-1):",
                  ss$num_ph1_subj_p))
    message(paste("Number of subjects (vaccine group, phase-2):",
                  ss$num_ph2_subj_v))
    message(paste("Number of subjects (placebo group, phase-2):",
                  ss$num_ph2_subj_p))

    # Number of events
    message(paste("Number of events (vaccine group, phase-1):",
                  ss$num_ph1_events_v))
    message(paste("Number of events (placebo group, phase-1):",
                  ss$num_ph1_events_p))
    message(paste("Number of events (vaccine group, phase-2):",
                  ss$num_ph2_events_v))
    message(paste("Number of events (placebo group, phase-2):",
                  ss$num_ph2_events_p))

    # Event proportions
    message(paste("Proportion of subjects with an event (vaccine group, phase-",
                  "1):", ss$prop_ph1_events_v))
    message(paste("Proportion of subjects with an event (placebo group, phase-",
                  "1):", ss$prop_ph1_events_p))
    message(paste("Proportion of subjects with an event (vaccine group, phase-",
                  "2):", ss$prop_ph2_events_v))
    message(paste("Proportion of subjects with an event (placebo group, phase-",
                  "2):", ss$prop_ph2_events_p))

  }

  invisible(ss)

}
