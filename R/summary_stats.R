#' Calculate summary statistics
#'
#' @description TO DO
#' @param dat A data object returned by `load_data`.
#' @param quietly Boolean. If true, output will not be printed.
#' @return A list containing the following: \itemize{
#'     \item{\code{one}: asdf}
#'     \item{\code{two}: asdf}
#'     \item{\code{three}: asdf}
#' }
#' @examples
#' print("to do")
#' @export
summary_stats <- function(dat, quietly=F) {

  # !!!!! Build this out

  if (!methods::is(dat,"dat_vaccine")) {
    stop(paste0("`dat` must be an object of class 'dat_vaccine' returned by lo",
                "ad_data()."))
  }

  ss <- list(
    "num_ph1_events_v" = sum(dat$v$delta==1),
    "num_ph1_events_p" = sum(dat$p$delta==1),
    "num_ph2_events_v" = sum(dat$v$z==1 & dat$v$delta==1),
    "num_ph2_events_p" = sum(dat$p$z==1 & dat$p$delta==1)
  )

  if (!quietly) {
    message(paste("Number of events (vaccine group, phase-1 subcohort):",
                  ss$num_ph1_events_v,"\n"))
    message(paste("Number of events (placebo group, phase-1 subcohort):",
                  ss$num_ph1_events_p,"\n"))
    message(paste("Number of events (vaccine group, phase-2 subcohort):",
                  ss$num_ph2_events_v,"\n"))
    message(paste("Number of events (placebo group, phase-2 subcohort):",
                  ss$num_ph2_events_p,"\n"))
  }

  return(ss) # !!!!! need to return "quietly" (i.e. don't print); use `invisible`?

}
