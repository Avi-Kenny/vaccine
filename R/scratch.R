
#
if (F) {

  # Generics
  fn2 <- function(y) UseMethod("fn2")
  fn2 <- function(y=44) {
    print(y)
  }
  fn2(22)
  fn2()

  # est_np() interface
  if (F) {
    dat <- load_data(999)
    ests <- est_np(dat=dat, cve=T, cr=T, ci_type="logit", cf_folds=0,
                   edge_corr="none")
  }

  # Misc
  dat <- list(
    a = c(1:10),
    # s = runif(5)
    s = c(runif(5),NA)
  )
  fn <- function(dat, x_out=seq(from=min(dat$s), to=max(dat$s), l=5)) {
    if (missing(x_out) && any(is.na(dat$s))) {
      x_out <- seq(from=min(dat$s, na.rm=T), to=max(dat$s, na.rm=T), l=101)
    }
    return(x_out)
  }
  fn(dat)
  fn(dat, x_out=c(88,99))

}
