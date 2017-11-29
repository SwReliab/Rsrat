#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

//' @rdname em
// [[Rcpp::export]]

List em_exp_emstep(NumericVector params, List data) {
  const double omega = params[0];
  const double rate = params[1];
  const int dsize = data["len"];
  NumericVector time = as<NumericVector>(data["time"]);
  NumericVector num = as<NumericVector>(data["fault"]);
  IntegerVector type = as<IntegerVector>(data["type"]);

  if (dsize != time.length() || dsize != num.length() || dsize != type.length()) {
    stop("Invalid data.");
  }

  double en1 = 0.0;
  double en2 = 0.0;
  double llf = 0.0;

  double t0 = 0.0;
  double t1 = time[0];
  double x1 = num[0];

  double tmp1, tmp2;

  // E-step
  if (x1 != 0.0) {
    tmp1 = 1 - exp(-rate*t1);
    tmp2 = 1.0/rate - (t1 + 1/rate) * exp(-rate*t1);
    en1 = x1;
    en2 = x1 * tmp2 / tmp1;
    llf = x1 * log(tmp1) - lgamma(x1+1.0);
  }
  if (type[0] == 1) {
    en1 += 1.0;
    en2 += t1;
    llf += log(rate) - rate*t1;
  }
  for (int j=1; j<dsize; j++) {
    if (time[j] != 0.0) {
      t0 = t1;
      t1 = t0 + time[j];
    }
    x1 = num[j];
    if (x1 != 0.0) {
      tmp1 = exp(-rate*t0) - exp(-rate*t1);
      tmp2 = (t0 + 1.0/rate) * exp(-rate*t0) - (t1 + 1.0/rate) * exp(-rate*t1);
      en1 += x1;
      en2 += x1 * tmp2 / tmp1;
      llf += x1 * log(tmp1) - lgamma(x1+1.0);
    }
    if (type[j] == 1) {
      en1 += 1.0;
      en2 += t1;
      llf += log(rate) - rate*t1;
    }
  }
  llf += log(omega) * en1;  // en1 is total number of faults
  en1 += omega * exp(-rate*t1);  // t1 is the last time
  en2 += omega * (t1 + 1.0/rate) * exp(-rate*t1); // t1 is the last time
  llf += - omega * (1.0 - exp(-rate*t1));

  // M-step
  double total = en1;
  double new_omega = en1;
  double new_rate = en1 / en2;

  return List::create(
    Named("param") = NumericVector::create(new_omega, new_rate),
    Named("pdiff") = NumericVector::create(new_omega - omega, new_rate - rate),
    Named("llf") = llf,
    Named("total") = total
  );
}
