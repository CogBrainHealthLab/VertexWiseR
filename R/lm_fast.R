#' Multiple-outcome linear regression t statistics
#'
#' Fits ordinary least squares models with one common design matrix using QR.
#' @param X Numeric design matrix with observations in rows and predictors in
#'   columns. Include the intercept explicitly, for example with model.matrix().
#' @param Y Numeric outcome matrix with observations in rows and outcomes in
#'   columns, or a numeric vector for a single outcome. Rows must align with X.
#' @param block_size Positive integer specifying outcomes per processing block.
#' @param tol Reciprocal-condition-number cutoff for the column-scaled design.
#'   This is not the same rank tolerance as in .lm.fit().
#' @param return_coefficients Logical. If TRUE, also return fitted coefficients
#'   in the original units of X. Defaults to FALSE.
#' @return By default, a numeric matrix of conventional coefficient t statistics
#'   for a null value of zero. If return_coefficients = TRUE, a list with matrices
#'   named t and coefficients. All matrices have predictors in rows and outcomes
#'   in columns, with names inherited from X and Y.
#' @details Requires finite data, full column rank and nrow(X) > ncol(X).
#'   Residual degrees of freedom are nrow(X) - ncol(X). No missing observations
#'   are removed. These are ordinary, not robust or mixed-model, standard errors.
#'   Exact or nearly perfect fits produce undefined or unstable t statistics.
#'   The QR decomposition is shared across outcomes within each call.
#' @importFrom Rcpp evalCpp
#' @export
#' @examples
#' set.seed(1)
#' X <- cbind(Intercept = 1, age = rnorm(100))
#' Y <- matrix(rnorm(100 * 5), nrow = 100)
#' t_values <- lm_fast(X, Y)
#' t_values["age", ]
#' fit <- lm_fast(X, Y, return_coefficients = TRUE)
#' fit$coefficients["age", ]
lm_fast <- function(X, Y, block_size = 2048L, tol = 1e-12,
                    return_coefficients = FALSE) {
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  if (!is.numeric(X) || !is.numeric(Y) || is.complex(X) || is.complex(Y)) {
    stop("X and Y must contain real numeric values.", call. = FALSE)
  }
  storage.mode(X) <- "double"
  storage.mode(Y) <- "double"
  if (ncol(X) < 1L || nrow(X) <= ncol(X) ||
      nrow(Y) != nrow(X) || ncol(Y) < 1L) {
    stop("Require nrow(X) == nrow(Y), nrow(X) > ncol(X), and nonempty columns.",
         call. = FALSE)
  }
  if (!is.numeric(block_size) || is.complex(block_size) ||
      length(block_size) != 1L || !is.finite(block_size) ||
      block_size < 1 || block_size > .Machine$integer.max ||
      block_size != floor(block_size)) {
    stop("block_size must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(tol) || is.complex(tol) || length(tol) != 1L ||
      !is.finite(tol) || tol <= 0 || tol >= 1) {
    stop("Require 0 < tol < 1.", call. = FALSE)
  }
  if (!is.logical(return_coefficients) || length(return_coefficients) != 1L ||
      is.na(return_coefficients)) {
    stop("return_coefficients must be TRUE or FALSE.", call. = FALSE)
  }
  # lm_multi_t_cpp is generated in R/RcppExports.R by compileAttributes().
  ans <- lm_multi_t_cpp(X, Y, as.integer(block_size), as.double(tol),
                        return_coefficients)
  dn <- list(colnames(X), colnames(Y))
  if (return_coefficients) {
    dimnames(ans$t) <- dn
    dimnames(ans$coefficients) <- dn
  } else {
    dimnames(ans) <- dn
  }
  ans
}
