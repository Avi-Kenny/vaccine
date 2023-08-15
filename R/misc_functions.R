#' Various elementary functions
#'
#' @param x Numeric input
#' @return Numeric output
#' @noRd
logit <- function(x) { log(x/(1-x)) }
expit <- function(x) { 1 / (1+exp(-x)) }
expit2 <- function(x) { (expit(x)-0.001)/0.998 }
logit2 <- function(x) { logit(0.001+0.998*x) }
log2 <- function(x) { log(x+0.001) }
exp2 <- function(x) { exp(x) - 0.001 }
deriv_expit <- function(x) { exp(x) / ((1+exp(x))^2) }
deriv_logit <- function(x) { 1 / (x-x^2) }
deriv_logit2 <- function(x) { 0.998*deriv_logit(0.001+0.998*x) }
deriv_log2 <- function(x) { 1 / (x+0.001) }



#' Alias for indicator (as.integer) function
#'
#' @noRd
In <- as.integer



#' Helper function for debugging; prints timestamps
#'
#' @noRd
chk <- function(num, msg="") {
  if (msg=="") {
    str <- paste0("Check ", num, ": ", Sys.time())
  } else {
    str <- paste0("Check ", num, " (", msg, "): ", Sys.time())
  }
  message(str)
}



#' Memoise a function
#'
#' @param fnc A function to be memoised
#' @return Memoised version of function
#' @note
#'   - This is a lightweight/faster version of the `memoise` function
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
#' @param dat_orig Dataset returned by `load_data`
#' @param grid_size A list, as specified in `est_np`
#' @return Dataset with values rounded
#' @noRd
create_grid <- function(dat_orig, grid_size, t_0) {

  d <- dat_orig
  grid <- list()
  grid$y <- round(seq(from=0, to=t_0, length.out=grid_size$y), 5)
  grid$y_ext <- round(seq(from=0, to=max(dat_orig$y),
                          by=(t_0/(grid_size$y-1))), 5)
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
#' @param dat_orig Dataset returned by `load_data`
#' @param grid A grid, returned by `create_grid`
#' @param grid_size A list, as specified in `est_np`
#' @return Dataset with values rounded
#' @noRd
round_dat <- function(dat_orig, grid, grid_size) {

  d <- dat_orig

  # Round `y` and `s`
  d$y <- sapply(d$y, function(y) { grid$y_ext[which.min(abs(grid$y_ext-y))] })
  d$s <- sapply(d$s, function(s) {
    if (is.na(s)) {
      return(NA)
    } else {
      return(grid$s[which.min(abs(grid$s-s))])
    }
  })

  # Round `x`
  for (i in c(1:length(d$x))) {
    x_col <- d$x[,i]
    if (length(unique(x_col))>grid_size$x) {
      d$x[,i] <- sapply(x_col, function(x) {
        grid$x[[i]][which.min(abs(grid$x[[i]]-x))]
      })
    }
  }

  # Not rounding weights for now

  return(d)

}



#' Subset dat_orig according to indices
#'
#' @param dat_orig Dataset returned by `load_data`
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
  )
  if (!is.null(dat_orig$spl)) {
    dat$spl <- dat_orig$spl[i,, drop=F]
  }
  if (!is.null(dat_orig$strata)) {
    dat$strata <- dat_orig$strata[i]
  }
  attr(dat, "n_orig") <- attr(dat_orig, "n_orig")
  attr(dat, "dim_x") <- attr(dat_orig, "dim_x")

  return(dat)

}



#' Convert dat_orig or dat to a data frame
#'
#' @param d Either dat_orig or dat
#' @return Data frame version of data object
#' @noRd
as_df <- function(d, strata=F) {
  df <- cbind(d$x, s=d$s, y=d$y, delta=d$delta, z=d$z, weights=d$weights)
  if (strata) { df$strata <- d$strata }
  return(df)
}



#' Copy of apply function
#'
#' @noRd
apply2 <- function (X, MARGIN, FUN, ..., simplify=TRUE) {
  FUN <- match.fun(FUN)
  simplify <- isTRUE(simplify)
  dl <- length(dim(X))
  if (!dl)
    stop("dim(X) must have a positive length")
  if (is.object(X))
    X <- if (dl == 2L)
      as.matrix(X)
  else as.array(X)
  d <- dim(X)
  dn <- dimnames(X)
  ds <- seq_len(dl)
  if (is.character(MARGIN)) {
    if (is.null(dnn <- names(dn)))
      stop("'X' must have named dimnames")
    MARGIN <- match(MARGIN, dnn)
    if (anyNA(MARGIN))
      stop("not all elements of 'MARGIN' are names of dimensions")
  }
  d.call <- d[-MARGIN]
  d.ans <- d[MARGIN]
  if (anyNA(d.call) || anyNA(d.ans))
    stop("'MARGIN' does not match dim(X)")
  s.call <- ds[-MARGIN]
  s.ans <- ds[MARGIN]
  dn.call <- dn[-MARGIN]
  dn.ans <- dn[MARGIN]
  d2 <- prod(d.ans)
  if (d2 == 0L) {
    newX <- array(vector(typeof(X), 1L), dim = c(prod(d.call),
                                                 1L))
    ans <- forceAndCall(1, FUN, if (length(d.call) < 2L) newX[,
                                                              1] else array(newX[, 1L], d.call, dn.call), ...)
    return(if (is.null(ans)) ans else if (length(d.ans) <
                                          2L) ans[1L][-1L] else array(ans, d.ans, dn.ans))
  }
  newX <- aperm(X, c(s.call, s.ans))
  dim(newX) <- c(prod(d.call), d2)
  ans <- vector("list", d2)
  if (length(d.call) < 2L) {
    if (length(dn.call))
      dimnames(newX) <- c(dn.call, list(NULL))
    for (i in 1L:d2) {
      tmp <- forceAndCall(1, FUN, newX[, i], ...)
      if (!is.null(tmp))
        ans[[i]] <- tmp
    }
  }
  else for (i in 1L:d2) {
    tmp <- forceAndCall(1, FUN, array(newX[, i], d.call,
                                      dn.call), ...)
    if (!is.null(tmp))
      ans[[i]] <- tmp
  }
  ans.list <- !simplify || is.recursive(ans[[1L]])
  l.ans <- length(ans[[1L]])
  ans.names <- names(ans[[1L]])
  if (!ans.list)
    ans.list <- any(lengths(ans) != l.ans)
  if (!ans.list && length(ans.names)) {
    all.same <- vapply(ans, function(x) identical(names(x),
                                                  ans.names), NA)
    if (!all(all.same))
      ans.names <- NULL
  }
  len.a <- if (ans.list)
    d2
  else length(ans <- unlist(ans, recursive = FALSE))
  if (length(MARGIN) == 1L && len.a == d2) {
    names(ans) <- if (length(dn.ans[[1L]]))
      dn.ans[[1L]]
    ans
  }
  else if (len.a == d2)
    array(ans, d.ans, dn.ans)
  else if (len.a && len.a%%d2 == 0L) {
    if (is.null(dn.ans))
      dn.ans <- vector(mode = "list", length(d.ans))
    dn1 <- list(ans.names)
    if (length(dn.call) && !is.null(n1 <- names(dn <- dn.call[1])) &&
        nzchar(n1) && length(ans.names) == length(dn[[1]]))
      names(dn1) <- n1
    dn.ans <- c(dn1, dn.ans)
    array(ans, c(len.a%/%d2, d.ans), if (!is.null(names(dn.ans)) ||
                                         !all(vapply(dn.ans, is.null, NA)))
      dn.ans)
  }
  else ans
}



#' Monotonize confidence limits (NP estimator)
#'
#' @param ci_lo A vector of confidence interval lower limits
#' @param ci_up A vector of confidence interval upper limits
#' @param dir Direction of monotonicity; one of c("decr", "incr")
#' @param type One of c("regular", "conservative")
#' @return A list of the form list(ci_lo=c(), ci_up=c())
#' @noRd
monotonize_cis <- function(ci_lo, ci_up, dir, type="regular") {

  # !!!!!
  if (F) {

    # # Set 1
    # dir <- "incr"
    # # type <- "conservative"
    # type <- "regular"
    # ci_lo <- c(1.1, 1.3, 1.5, 1.2, 1.6, 1.8, 1.7, 2.0)
    # ci_up <- c(2.1, 2.5, 2.4, 2.2, 2.6, 2.8, 2.9, 2.5)

    # Set 2
    dir <- "decr"
    # type <- "conservative"
    type <- "regular"
    ci_lo <- rev(c(1.1, 1.3, 1.5, 1.2, 1.6, 1.8, 1.7, 2.0))
    ci_up <- rev(c(2.1, 2.5, 2.4, 2.2, 2.6, 2.8, 2.9, 2.5))

    ci_lo_orig <- ci_lo; ci_up_orig <- ci_up

    df_plot <- data.frame(
      x = rep(seq(1:8), 4),
      y = c(ci_lo, ci_up, ci_lo_orig, ci_up_orig),
      which = rep(c("lower", "upper", "lower", "upper"), each=8),
      type = rep(c("updated", "old"), each=16)
    )
    ggplot(df_plot, aes(x=x, y=y, color=factor(which), linetype=type)) +
      geom_line() +
      labs(title=paste0(dir, "; ", type))

  }

  # This helper function returns the "least nondecreasing majorant"
  lnm <- function(y) {
    val <- y[1]
    for (i in c(2:length(y))) {
      if (!is.na(y[i]) && !is.na(val) && y[i]<val) { y[i] <- val }
      val <- y[i]
    }
    return(y)
  }

  if (type=="regular") {
    if (dir=="decr") {
      ci_lo <- rev(lnm(rev(ci_lo)))
      ci_up <- -1*lnm(-1*ci_up)
    } else {
      ci_lo <- lnm(ci_lo)
      ci_up <- -1*rev(lnm(rev(-1*ci_up)))
    }
  } else if (type=="conservative") {
    if (dir=="decr") {
      ci_lo <- -1*lnm(-1*ci_lo)
      ci_up <- rev(lnm(rev(ci_up)))
    } else {
      ci_lo <- -1*rev(lnm(rev(-1*ci_lo)))
      ci_up <- lnm(ci_up)
    }
  }

  return(list(ci_lo=ci_lo, ci_up=ci_up))

}
