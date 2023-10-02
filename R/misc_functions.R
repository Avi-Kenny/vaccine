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
#' @param dat Dataset returned by `load_data` (possibly filtered)
#' @param grid_size A list, as specified in `params_ce_np`
#' @return Dataset with values rounded
#' @noRd
create_grid <- function(dat, grid_size, t_0) {

  grid <- list()
  grid$y <- round(seq(from=0, to=t_0, length.out=grid_size$y), 5)
  grid$y_ext <- round(seq(from=0, to=max(dat$y), by=(t_0/(grid_size$y-1))), 5)
  if (!(t_0 %in% grid$y)) { grid$y <- sort(c(grid$y, t_0)) }
  grid$s <- round(seq(from=0, to=1, length.out=grid_size$s), 5)
  grid$x <- lapply(c(1:attr(dat,"dim_x")), function(i) {
    x_col <- as.numeric(dat[,i])
    if (length(unique(x_col[!is.na(x_col)]))>grid_size$x) {
      return(round(seq(from=min(x_col, na.rm=T), to=max(x_col, na.rm=T),
                       length.out=grid_size$x), 5))
    } else {
      return(NA)
    }
  })

  return(grid)

}



#' Round data values
#'
#' @param dat Dataset returned by `load_data` (possibly filtered)
#' @param grid A grid, returned by `create_grid`
#' @param grid_size A list, as specified in `est_np`
#' @return Dataset with values rounded
#' @noRd
round_dat <- function(dat, grid, grid_size) {

  # Round `y` and `s`
  dat$y <- sapply(dat$y, function(y) {
    grid$y_ext[which.min(abs(grid$y_ext-y))]
  })
  dat$s <- sapply(dat$s, function(s) {
    if (is.na(s)) {
      return(NA)
    } else {
      return(grid$s[which.min(abs(grid$s-s))])
    }
  })

  # Round `x`
  for (i in c(1:attr(dat,"dim_x"))) {
    x_col <- as.numeric(dat[,i])
    if (length(unique(x_col[!is.na(x_col)]))>grid_size$x) {
      dat[,i] <- sapply(x_col, function(x) {
        if (is.na(x)) {
          return(NA)
        } else {
          return(grid$x[[i]][which.min(abs(grid$x[[i]]-x))])
        }
      })
    }
  }

  # !!!!! Note: not rounding weights (for now)

  return(dat)

}



#' Find dataframe index of a row
#'
#' @param vec A numeric vector
#' @param df A dataframe to search
#' @return The index of the row of `df` at which `vec` is found
#' @noRd
find_index <- function(vec, df) {
  vec <- as.numeric(vec)
  r <- list()
  for (i in c(1:length(vec))) { r[[i]] <- which(abs(vec[i]-df[,i])<1e-8) }
  index <- Reduce(intersect, r)
  if (length(index)!=1) {
    if (length(index)==0) { stop("length(index)==0") }
    if (length(index)>1) { stop("length(index)>1") }
  }
  return(index)
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
