//
#include <RcppArmadillo.h>
#include <math.h>

//#include <RcppEigen.h>

using namespace arma;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

double shooting(double z, double r)
{
  if (z > r) {
    return z - r;
  } else if (z < -r) {
    return z + r;
  } else {
    return .0;
  }
}

// [[Rcpp::export]]
Rcpp::List gaussian_elastic_arma(arma::vec &y, arma::mat &x, arma::vec &lambdaVec, double &alpha,
                                 arma::uword &maxit, double &tol)
{
  arma::uword n = x.n_rows;
  double n_double = (double)n;
  arma::uword d = x.n_cols;
  double d_double = (double)d;
  double b0 = 0;
  arma::vec b(d, fill::zeros);
  arma::uword lambdaLength = lambdaVec.n_elem;
  //
  arma::vec b0Vec(lambdaLength, fill::zeros);
  arma::mat bMat(d, lambdaLength, fill::zeros);
  //
  double sy = arma::accu(y);
  arma::vec yj(n, fill::zeros);
  double stopCrit;
  for (arma::uword lm = 0; lm < lambdaLength; lm++) {
    double lambda = lambdaVec(lm);
    for (arma::uword i = 0; i < maxit; i++) {
      b0 = (sy - arma::accu(sum(x*b)))/n_double;
      for (arma::uword j = 0; j < d; j++) {
        yj = b0 + x*b - x.col(j)*b(j);
        b(j) = shooting(arma::as_scalar((x.col(j).t()*(y - yj)))/n_double, lambda*alpha)/(1. + lambda*(1. - alpha));
      }
      // Compute the mean absolute difference
      stopCrit = (std::abs(b0Vec(lm) - b0) + arma::accu(arma::abs(bMat.col(lm) - b)))/(d_double + 1);
      // Update coefficients
      b0Vec(lm) = b0;
      bMat.col(lm) = b;
      // Check stopping criterion
      if (stopCrit < tol) { i += maxit; }
    }
  }
  return List::create(
    Named("b0") = wrap(b0),
    Named("b") = wrap(bMat),
    Named("lambda") = wrap(lambdaVec)
  );
}

