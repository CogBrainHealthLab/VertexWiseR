#' Threshold-Free Cluster Enhancement (TFCE) in 3D
#'
#' @param data surface data
#' @param tail `1`,`-1` or `2` for positive, negative or two-tailed
#' @return TFCE object
#' @useDynLib VertexWiseR, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @export
TFCE <- function(data, tail, edgelist) {
  data <- as.numeric(data)
  tail <- as.integer(tail)
  edgelist <- as.matrix(edgelist)
  storage.mode(edgelist) <- "integer"
  return(TFCE_cpp_impl(data, tail, edgelist))
}