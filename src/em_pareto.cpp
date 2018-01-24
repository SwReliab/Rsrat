#include <Rcpp.h>
#include <cmath>

#include "numlib.h"
#include "inv_log_psi.h"

using namespace Rcpp;

//' @rdname em
// [[Rcpp::export]]

List em_pareto_emstep(NumericVector params, List data) {
  const double omega = params[0];
  const double shape = params[1];
  const double scale = params[2];
  const int dsize = data["len"];
  NumericVector time = as<NumericVector>(data["time"]);
  NumericVector num = as<NumericVector>(data["fault"]);
  IntegerVector type = as<IntegerVector>(data["type"]);

  if (dsize != time.length() || dsize != num.length() || dsize != type.length()) {
    stop("Invalid data.");
  }

  double t = time[0];
  double x = num[0];

  double en1 = 0.0;
  double en2 = 0.0;
  double en3 = 0.0;
  double llf = 0.0;

  double g00 = 1.0;
  double g01 = shape/scale;
  double g02 = psi(shape) - log(scale);
  double g10 = pow(scale/(scale+t), shape); // 1.0 - pareto_cdf(t, shape, scale);
  double g11 = shape / (scale + t) * g10;
  double g12 = (psi(shape) - log(scale + t)) * g10;

  if (x != 0.0) {
    double tmp1 = g00 - g10;
    double tmp2 = g01 - g11;
    double tmp3 = g02 - g12;
    en1 += x;
    en2 += x * tmp2 / tmp1;
    en3 += x * tmp3 / tmp1;
    llf += x * log(tmp1) - lgamma(x+1);
  }
  if (type[0] == 1) {
    en1 += 1.0;
    en2 += (shape+1.0)/(scale+t);
    en3 += psi(shape+1.0) - log(scale + t);
    llf += log(shape) + shape * log(scale) - (shape+1) * log(scale + t); // log(pareto_pdf(t, shape, scale));
  }
  for (int j=1; j<dsize; j++) {
    x = num[j];
    if (time[j] != 0.0) {
      t += time[j];
      g00 = g10;
      g01 = g11;
      g02 = g12;
      g10 = pow(scale/(scale+t), shape); // 1.0 - pareto_cdf(t, shape, scale);
      g11 = shape / (scale + t) * g10;
      g12 = (psi(shape) - log(scale + t)) * g10;
    }
    if (x != 0.0) {
      double tmp1 = g00 - g10;
      double tmp2 = g01 - g11;
      double tmp3 = g02 - g12;
      en1 += x;
      en2 += x * tmp2 / tmp1;
      en3 += x * tmp3 / tmp1;
      llf += x * log(tmp1) - lgamma(x+1);
    }
    if (type[j] == 1) {
      en1 += 1.0;
      en2 += (shape+1.0)/(scale+t);
      en3 += psi(shape+1.0) - log(scale + t);
      llf += log(shape) + shape * log(scale) - (shape+1) * log(scale + t); // log(pareto_pdf(t, shape, scale));
    }
  }
  llf += log(omega) * en1;  // en1 is total number of faults
  en1 += omega * g10;  // g10 is the last time
  en2 += omega * g11;  // g11 is the last time
  en3 += omega * g12;  // g12 is the last time
  llf += - omega * (1.0 - g10);

  double total = en1;
  double new_omega =  en1;
  double new_shape = inv_log_psi(log(en2/en1)-en3/en1);
  double new_scale = new_shape * en1 / en2;

  return List::create(
    Named("param") = NumericVector::create(new_omega, new_shape, new_scale),
    Named("pdiff") = NumericVector::create(new_omega - omega, new_shape - shape, new_scale - scale),
    Named("llf") = llf,
    Named("total") = total
  );
}
