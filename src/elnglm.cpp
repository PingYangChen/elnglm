#include "header.h"

double shooting(const double &z, const double &r)
{
  if (z > r) {
    return z - r;
  } else if (z < -r) {
    return z + r;
  } else {
    return .0;
  }
}

arma::vec prob_value(const arma::mat &x, const double &b0, const arma::vec &b) 
{
  arma::vec val = 1./(1. + arma::exp(-(b0 + x*b)));
  arma::vec onev(x.n_rows, fill::ones);
  val = arma::max(val, onev*1e-5);
  val = arma::min(val, onev*(1.- 1e-5));
  return val;
}

void calcWZ(arma::vec &w, arma::vec &z, 
            const arma::vec &y, const arma::mat &x, const double &b0, const arma::vec &b)
{
  arma::vec pv = prob_value(x, b0, b);
  //w = pv % (1. - pv);
  w.fill(.25);
  z = b0 + x*b + (y - pv)/w;
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
  arma::vec r(n, fill::zeros);
  arma::vec tmp(n, fill::zeros);
  double lambda;
  double stopCrit;
  for (arma::uword lm = 0; lm < lambdaLength; lm++) {
    lambda = lambdaVec(lm);
    for (arma::uword i = 0; i < maxit; i++) {
      r = y - b0 - x*b;
      b0 = arma::accu(r + b0)/n_double;
      for (arma::uword j = 0; j < d; j++) {
        tmp = x.col(j)*b(j);
        b(j) = shooting(arma::as_scalar(x.col(j).t()*(r + tmp))/n_double, lambda*alpha)/(1. + lambda*(1. - alpha));
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
    Named("b0") = wrap(b0Vec),
    Named("b") = wrap(bMat),
    Named("lambda") = wrap(lambdaVec)
  );
}

// [[Rcpp::export]]
Rcpp::List binomial_elastic_arma(arma::vec &y, arma::mat &x, arma::vec &lambdaVec, double &alpha,
                                 arma::uword &maxit, double &tol)
{
  arma::uword n = x.n_rows;
  //double n_double = (double)n;
  arma::uword d = x.n_cols;
  double d_double = (double)d;
  arma::vec w(n);
  arma::vec z(n);
  double b0 = 0.;
  arma::vec b(d, fill::zeros);
  arma::uword lambdaLength = lambdaVec.n_elem;
  //
  arma::vec b0Vec(lambdaLength, fill::zeros);
  arma::mat bMat(d, lambdaLength, fill::zeros);
  //
  arma::vec r(n, fill::zeros);
  arma::vec tmp(n, fill::zeros);
  double lambda;
  double stopCrit;
  for (arma::uword lm = 0; lm < lambdaLength; lm++) {
    lambda = lambdaVec(lm);
    for (arma::uword i = 0; i < maxit; i++) {
      calcWZ(w, z, y, x, b0, b);
      r = z - b0 - x*b;
      b0 = arma::accu(w % (r + b0))/arma::accu(w);
      for (arma::uword j = 0; j < d; j++) {
        tmp = x.col(j)*b(j);
        b(j) = shooting(arma::as_scalar((w % x.col(j)).t()*(r + tmp)), lambda*alpha)/(arma::accu(w % x.col(j) % x.col(j)) + lambda*(1 - alpha));
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
    Named("b0") = wrap(b0Vec),
    Named("b") = wrap(bMat),
    Named("lambda") = wrap(lambdaVec)
  );
}

// [[Rcpp::export]]
Rcpp::List multinomial_elastic_arma(arma::mat &y, arma::mat &x, arma::vec &lambdaVec, double &alpha,
                                    arma::uword &maxit, double &tol)
{
  arma::uword n = x.n_rows;
  //double n_double = (double)n;
  arma::uword d = x.n_cols;
  double d_double = (double)d;
  arma::uword yd = y.n_cols;
  arma::vec w(n);
  arma::vec z(n);
  arma::vec b0(yd, fill::zeros);
  arma::mat b(d, yd, fill::zeros);
  arma::uword lambdaLength = lambdaVec.n_elem;
  //
  arma::mat b0Mat(yd, lambdaLength, fill::zeros);
  arma::cube bCube(d, yd, lambdaLength, fill::zeros);
  //
  arma::vec r(n, fill::zeros);
  arma::vec tmp(n, fill::zeros);
  double lambda;
  double stopCrit;
  for (arma::uword lm = 0; lm < lambdaLength; lm++) {
    lambda = lambdaVec(lm);
    for (arma::uword m = 0; m < yd; m++) {
      for (arma::uword i = 0; i < maxit; i++) {
        calcWZ(w, z, y.col(m), x, b0(m), b.col(m));
        r = z - b0(m) - x*b.col(m);
        b0(m) = arma::accu(w % (r + b0(m)))/arma::accu(w);
        for (arma::uword j = 0; j < d; j++) {
          tmp = x.col(j)*b(j,m);
          b(j,m) = shooting(arma::as_scalar((w % x.col(j)).t()*(r + tmp)), lambda*alpha)/(arma::accu(w % x.col(j) % x.col(j)) + lambda*(1 - alpha));
        }
        // Compute the mean absolute difference
        stopCrit = (std::abs(b0Mat(m,lm) - b0(m)) + arma::accu(arma::abs(bCube.slice(lm).col(m) - b.col(m))))/(d_double + 1);
        // Update coefficients
        b0Mat(m,lm) = b0(m);
        bCube.slice(lm).col(m) = b.col(m);
        //bCube.subcube(0,m,lm,d-1,m,lm) = b.col(m);
        // Check stopping criterion
        if (stopCrit < tol) { i += maxit; }
      }
    }
  }
  return List::create(
    Named("b0") = wrap(b0Mat),
    Named("b") = wrap(bCube),
    Named("lambda") = wrap(lambdaVec)
  );
}
