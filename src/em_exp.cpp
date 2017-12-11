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

namespace expsrm {
  double func_barFi(double t, double rate) {
    return exp(-rate*t);
  }

  double func_barH1i(double t, double rate) {
    return (t + 1/rate) * exp(-rate*t);
  }
}

//' @rdname em
// [[Rcpp::export]]

List em_exp_emstep2(NumericVector params, List data) {
  const double omega = params[0];
  const double rate = params[1];
  const int dsize = data["len"];
  NumericVector time = as<NumericVector>(data["time"]);
  NumericVector num = as<NumericVector>(data["fault"]);
  IntegerVector type = as<IntegerVector>(data["type"]);

  if (dsize != time.length() || dsize != num.length() || dsize != type.length()) {
    stop("Invalid data.");
  }

  double nn = 0;
  double en1 = 0;
  double en2 = 0;
  double llf = 0;

  double t = 0;
  double prev_barFi = 1;
  double prev_barH1i = 1/rate;

  for (int i=0; i<dsize; i++) {
    t += time[i];
    double barFi = expsrm::func_barFi(t, rate);
    double barH1i = expsrm::func_barH1i(t, rate);

    double x = num[i];
    if (x != 0) {
      double tmp1 = prev_barFi - barFi;
      double tmp2 = prev_barH1i - barH1i;
      nn += x;
      en1 += x;
      en2 += x * tmp2 / tmp1;
      llf += x * log(tmp1) - lgamma(x+1);
    }
    if (type[i] == 1) {
      nn += 1;
      en1 += 1;
      en2 += t;
      llf += R::dexp(t, rate, true);
    }
    prev_barFi = barFi;
    prev_barH1i = barH1i;
  }
  llf += nn * log(omega) - omega * (1 - prev_barFi);
  en1 += omega * prev_barFi;
  en2 += omega * prev_barH1i;

  double new_omega = en1;
  double new_rate = en1 / en2;
  double total = nn + omega * exp(-new_rate*t);

  return List::create(
    Named("param") = NumericVector::create(new_omega, new_rate),
    Named("pdiff") = NumericVector::create(new_omega - omega, new_rate - rate),
    Named("llf") = llf,
    Named("total") = total
  );
}
