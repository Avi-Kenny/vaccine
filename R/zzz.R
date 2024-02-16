
.onAttach <- function(libname, pkgname) {
  vrs <- utils::packageVersion("vaccine")
  packageStartupMessage(paste0(
    "vaccine (version ", vrs, ").\nType ?vaccine to get started."
  ))
}
