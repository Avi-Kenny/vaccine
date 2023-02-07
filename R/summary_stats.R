#' Calculate summary statistics
#'
#' @description TO DO
#' @param dat A data object returned by `load_data`
#' @return A list containing the following: \itemize{
#'     \item{\code{one}: asdf}
#'     \item{\code{two}: asdf}
#'     \item{\code{three}: asdf}
#' }
#' @examples
#' print("to do")
#' @export
summary_stats <- function(dat, print=F) {

  if (class(dat)!="dat_vaccine") {
    stop(paste0("`dat` must be an object of class 'dat_vaccine' returned by lo",
                "ad_data()."))
  }

  ss <- list(
    "num_ph1_events_vc" = sum(df_vc$delta==1),
    "num_ph1_events_pl" = sum(df_pl$delta==1),
    "num_ph2_events_vc" = sum(df_vc$z==1 & df_vc$delta==1),
    "num_ph2_events_pl" = sum(df_pl$z==1 & df_pl$delta==1)
  )

  if (print) {
    message(paste("Number of events (vaccine group, phase-1 subcohort):",
                  ss$num_ph1_events_vc,"\n"))
    message(paste("Number of events (placebo group, phase-1 subcohort):",
                  ss$num_ph1_events_pl,"\n"))
    message(paste("Number of events (vaccine group, phase-2 subcohort):",
                  ss$num_ph2_events_vc,"\n"))
    message(paste("Number of events (placebo group, phase-2 subcohort):",
                  ss$num_ph2_events_pl,"\n"))
  }

  return(ss) # !!!!! need to return "quietly" (i.e. don't print)

}
