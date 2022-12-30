#' Various elementary functions
#'
#' @param x Numeric input
#' @return Numeric output
#' @noRd
logit <- function(x) { log(x/(1-x)) }
expit <- function(x) { 1 / (1+exp(-x)) }
deriv_expit <- function(x) { exp(x) / ((1+exp(x))^2) }
deriv_logit <- function(x) { 1 / (x-x^2) }



#' Alias for indicator (as.integer) function
#'
#' @noRd
In <- as.integer



#' Memoise a function
#'
#' @param fnc A function to be memoised
#' @return Memoised version of function
#' @note
#'   - This is a lightweight/faster version of the memoise() function
#' @noRd
memoise2 <- function(fnc) {

  htab <- new.env()
  ..new_fnc <- function() {
    ..e <- parent.env(environment())
    ..mc <- lapply(as.list(match.call())[-1L], eval, parent.frame())
    key <- rlang::hash(..mc)
    val <- ..e$htab[[key]]
    if (is.null(val)) {
      val <- do.call(..e$fnc, ..mc)
      ..e$htab[[key]] <- val
    }
    return(val)
  }

  # Set formals and set up environment
  formals(..new_fnc) <- formals(fnc)
  f_env <- new.env(parent=environment(fnc))
  f_env$arg_names <- names(formals(fnc))
  f_env$htab <- htab
  f_env$fnc <- fnc
  environment(..new_fnc) <- f_env

  return(..new_fnc)

}



#' Round data values
#'
#' @param dat_orig Dataset returned by load_data()
#' @param grid_size A list, as specified in est_np()
#' @return Dataset with values rounded
#' @noRd
create_grid <- function(dat_orig, grid_size, t_0) {

  d <- dat_orig
  grid <- list()
  grid$y <- round(seq(from=0, to=t_0, length.out=grid_size$y), 5)
  if (!(t_0 %in% grid$y)) { grid$y <- sort(c(grid$y, t_0)) }
  grid$s <- round(seq(from=0, to=1, length.out=grid_size$s), 5)
  grid$x <- lapply(c(1:length(d$x)), function(i) {
    x_col <- d$x[,i]
    if (length(unique(x_col))>grid_size$x) {
      return(round(seq(from=min(x_col), to=max(x_col),
                       length.out=grid_size$x), 5))
    } else {
      return(NA)
    }
  })

  return(grid)

}



#' Round data values
#'
#' @param dat_orig Dataset returned by load_data()
#' @param grid A grid, returned by create_grid()
#' @param grid_size A list, as specified in est_np()
#' @return Dataset with values rounded
#' @noRd
round_dat <- function(dat_orig, grid, grid_size) {

  d <- dat_orig

  # Round `y` and `s`
  d$y <- sapply(d$y, function(y) { grid$y[which.min(abs(grid$y-y))] })
  d$s <- sapply(d$s, function(s) {
    if (is.na(s)) {
      return(NA)
    } else {
      return(grid$s[which.min(abs(grid$s-s))])
    }
  })

  # Round `x`
  for (i in c(1:length(d$x))) {
    x_col <- d$x[i,]
    if (length(unique(x_col))>grid_size$x) {
      d$x[i,] <- sapply(x_col, function(x) {
        grid$x[[i]][which.min(abs(grid$x[[i]]-x))]
      })
    }
  }

  # Not rounding weights for now

  return(d)

}



#' Subset dat_orig according to indices
#'
#' @param dat_orig Dataset returned by load_data()
#' @param indices Indices to filter dataset by
#' @return Filtered subsample of dataset
#' @noRd
ss <- function(dat_orig, indices) {

  i <- indices
  dat <- list(
    y = dat_orig$y[i],
    delta = dat_orig$delta[i],
    s = dat_orig$s[i],
    x = dat_orig$x[i,, drop=F],
    weights = dat_orig$weights[i],
    z = dat_orig$z[i]
    # strata = dat_orig$strata[i] # !!!!!
  )
  attr(dat, "n_orig") <- attr(dat_orig, "n_orig")
  attr(dat, "dim_x") <- attr(dat_orig, "dim_x")

  return(dat)

}



#' Convert dat_orig or dat to a data frame
#'
#' @param d Either dat_orig or dat
#' @return Data frame version of data object
#' @noRd
as_df <- function(d) {
  cbind(d$x, s=d$s, y=d$y, delta=d$delta, z=d$z, weights=d$weights)
}
