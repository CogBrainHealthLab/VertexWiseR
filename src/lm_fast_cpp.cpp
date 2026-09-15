// Place in src/cpp.cpp. Generate bindings with Rcpp::compileAttributes().
// Package builds require CXX_STD = CXX14 in Makevars and Makevars.win.
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp14)]]

// [[Rcpp::export]]
Rcpp::RObject lm_multi_t_cpp(const arma::mat& X, const arma::mat& Y,
                         int block_size = 2048, double tol = 1e-12,
                         bool return_coefficients = false) {
  const arma::uword n = X.n_rows, p = X.n_cols, m = Y.n_cols;
  if (p == 0 || n <= p || Y.n_rows != n || m == 0)
    Rcpp::stop("Require nrow(X) == nrow(Y), n > p, p > 0 and at least one outcome.");
  if (block_size < 1 || !std::isfinite(tol) || tol <= 0 || tol >= 1)
    Rcpp::stop("Require block_size >= 1 and 0 < tol < 1.");
  if (!X.is_finite() || !Y.is_finite())
    Rcpp::stop("X and Y must contain only finite values.");

  // Positive column scaling improves numerical conditioning and leaves t unchanged.
  arma::mat Z = X;
  arma::vec scales(p);
  for (arma::uword j = 0; j < p; ++j) {
    const double scale = arma::norm(Z.col(j), 2);
    if (!std::isfinite(scale) || scale == 0)
      Rcpp::stop("X contains a zero or excessively large column.");
    Z.col(j) /= scale;
    scales(j) = scale;
  }
  arma::mat Q, R;
  if (!arma::qr_econ(Q, R, Z)) Rcpp::stop("QR decomposition failed.");
  const double rc = arma::rcond(R);
  if (!std::isfinite(rc) || rc <= tol)
    Rcpp::stop("X is rank deficient or too ill-conditioned; inspect the design.");

  // diag((Z transpose Z)^-1) = row sums of squared inverse(R).
  const arma::mat I = arma::eye<arma::mat>(p, p);
  arma::mat Ri;
  if (!arma::solve(Ri, arma::trimatu(R), I, arma::solve_opts::no_approx))
    Rcpp::stop("Triangular solve failed.");
  const arma::vec se_factor = arma::sqrt(arma::sum(arma::square(Ri), 1));
  arma::mat result(p, m);
  arma::mat coefficients;
  if (return_coefficients) coefficients.set_size(p, m);
  const double df = static_cast<double>(n - p);
  const arma::uword bs = static_cast<arma::uword>(block_size);

  for (arma::uword start = 0; start < m; start += bs) {
    Rcpp::checkUserInterrupt();
    const arma::uword end = std::min(m - 1, start + bs - 1);
    const arma::mat rhs = Q.t() * Y.cols(start, end);
    arma::mat beta;
    if (!arma::solve(beta, arma::trimatu(R), rhs, arma::solve_opts::no_approx))
      Rcpp::stop("Coefficient solve failed.");
    // Explicit residuals avoid subtraction of two nearly equal sums of squares.
    const arma::mat residual = Y.cols(start, end) - Z * beta;
    const arma::rowvec sigma = arma::sqrt(arma::sum(arma::square(residual), 0) / df);
    if (return_coefficients) {
      // Undo design-column scaling to return coefficients in original X units.
      arma::mat original_beta = beta;
      original_beta.each_col() /= scales;
      coefficients.cols(start, end) = original_beta;
    }
    beta.each_col() /= se_factor;
    beta.each_row() /= sigma;
    result.cols(start, end) = beta;
  }
  if (return_coefficients) {
    return Rcpp::List::create(Rcpp::Named("t") = result,
                              Rcpp::Named("coefficients") = coefficients);
  }
  return Rcpp::wrap(result);
}
