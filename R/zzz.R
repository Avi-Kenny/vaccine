
.onAttach <- function(libname, pkgname) {
  vrs <- utils::packageVersion("vaccine")
  packageStartupMessage(paste0(
    "vaccine (version ", vrs, ").\nType ?vaccine to get started."
  ))
}

.vaccine_env <- new.env(parent=emptyenv())
.onLoad <- function(libname, pkgname) {
  .vaccine_env$return_IF_vec <- FALSE
}
