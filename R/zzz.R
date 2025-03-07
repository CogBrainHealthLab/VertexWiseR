#' @importFrom utils packageVersion

.onLoad <- function(libname, pkgname) {
  if (!requireNamespace("reticulate", quietly = TRUE) || gsub(pattern = "\\.", "", utils::packageVersion("reticulate")) > "1410") {
    warning("This package works with reticulate version 1.41.0 or lower. It was not tested for versions above.\n")
  }
}