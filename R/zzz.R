#' @importFrom utils packageVersion

.onLoad <- function(libname, pkgname) {
  #check for reticulate version
  rv=unlist(utils::packageVersion("reticulate"))
  #only for three levels
  rv=as.numeric(paste(rv[1:3], collapse = ""))
    if (gsub(rv, pattern = "\\.", "", ) > "1410") {
    warning("This package works with reticulate version 1.41.0 or lower. It was not tested for versions above.\n")
  }
}
