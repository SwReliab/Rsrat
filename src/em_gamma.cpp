#include <Rcpp.h>
#include <cmath>

#include "numlib.h"
#include "gauss_inte.h"
#include "inv_log_psi.h"

using namespace Rcpp;

#define MAX_NN 100
#define ZERO 1.0e-10

static int n;
static double x[MAX_NN];
static double w[MAX_NN];
static double fx[MAX_NN];
static double fv[MAX_NN];

template < typename T > std::string to_string(const T& n) {
  std::ostringstream stm ;
  stm << n ;
  return stm.str();
}

double em_gamma_int(double t0, double t1, double shape, double rate) {
  double c = gauss_inte_fx(n, x, t0, t1, fx);
  for (int i=0; i<n; i++) {
    double y = rate * fx[i];
    double tmp = log(rate) + (shape-1.0) * log(y) - y - lgamma(shape);
    fv[i] = log(fx[i]) * exp(tmp);
  }
  return gauss_inte_fv(n, w, c, fv);
}

//' @rdname em
// [[Rcpp::export]]

List em_gamma_emstep(NumericVector params, List data, int divide = 15, double eps = 1.0e-10) {
  const double omega = params[0];
  const double shape = params[1];
  const double rate = params[2];
  const int dsize = data["len"];
  NumericVector time = as<NumericVector>(data["time"]);
  NumericVector num = as<NumericVector>(data["fault"]);
  IntegerVector type = as<IntegerVector>(data["type"]);

  if (dsize != time.length() || dsize != num.length() || dsize != type.length()) {
    stop("Invalid data.");
  }

  if (divide > MAX_NN) {
    stop("divide should be lower than " + to_string(MAX_NN));
  }
  n = divide; // n is a global variable.
  gauss_inte_w(divide, x, w, eps); // x, w are global variables.

  // E-step
  double a0 = lgamma(shape);
  double a1 = lgamma(shape + 1.0);
  double en1 = 0.0;
  double en2 = 0.0;
  double en3 = 0.0;
  double llf = 0.0;

  double t0 = 0.0;
  double t1 = time[0];
  double x1 = num[0];
  double gam10 = 1.0;
  double gam11 = 1.0;
  double gam20 = q_gamma(shape, rate*t1, a0);
  double gam21 = q_gamma(shape+1.0, rate*t1, a1);

  double tmp1, tmp2, tmp3, tmp4;

  tmp3 = em_gamma_int(ZERO, t1, shape, rate);
  tmp4 = tmp3;
  if (x1 != 0.0) {
    tmp1 = gam10 - gam20;
    tmp2 = (shape/rate) * (gam11 - gam21);
    en1 += x1;
    en2 += x1 * tmp2 / tmp1;
    en3 += x1 * tmp3 / tmp1;
    llf += x1 * log(tmp1) - lgamma(x1+1.0);
  }
  if (type[0] == 1) {
    en1 += 1.0;
    en2 += t1;
    en3 += log(t1);
    llf += shape*log(rate) + (shape-1.0)*log(t1) - rate*t1 - lgamma(shape);
  }
  for (int j=1; j<dsize; j++) {
    x1 = num[j];
    if (time[j] != 0.0) {
      t0 = t1;
      t1 = t0 + time[j];
      gam10 = gam20;
      gam11 = gam21;
      gam20 = q_gamma(shape, rate*t1, a0);
      gam21 = q_gamma(shape+1, rate*t1, a1);
      tmp3 = em_gamma_int(t0, t1, shape, rate);
      tmp4 += tmp3;
    }
    if (x1 != 0.0) {
      tmp1 = gam10 - gam20;
      tmp2 = (shape/rate) * (gam11 - gam21);
      en1 += x1;
      en2 += x1 * tmp2 / tmp1;
      en3 += x1 * tmp3 / tmp1;
      llf += x1 * log(tmp1) - lgamma(x1+1.0);
    }
    if (type[j] == 1) {
      en1 += 1.0;
      en2 += t1;
      en3 += log(t1);
      llf += shape*log(rate) + (shape-1.0)*log(t1) - rate*t1 - lgamma(shape);
    }
  }
  llf += log(omega) * en1;  // en1 is total number of faults
  en1 += omega * gam20;  // gam20 is the last time
  en2 += omega * (shape/rate) * gam21;  // gam21 is the last time
  en3 += omega * (psi(shape) - log(rate) - tmp4);
  llf += - omega * (1.0 - gam20);

  // M-step
  double total = en1;
  double new_omega =  en1;
  double new_shape = inv_log_psi(log(en2/en1)-en3/en1);
  double new_rate = new_shape * en1 / en2;

  return List::create(
    Named("param") = NumericVector::create(new_omega, new_shape, new_rate),
    Named("pdiff") = NumericVector::create(new_omega - omega, new_shape - shape, new_rate - rate),
    Named("llf") = llf,
    Named("total") = total
  );

}
