#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

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
      llf += log(rate) - rate * t; // 'R::dexp(t, rate, true)' did not work correctly.
    }
    prev_barFi = barFi;
    prev_barH1i = barH1i;
  }
  llf += nn * log(omega) - omega * (1 - prev_barFi);
  en1 += omega * prev_barFi;
  en2 += omega * prev_barH1i;

  double new_omega = en1;
  double new_rate = en1 / en2;
  double total = en1; // nn + omega * exp(-new_rate*t);

  return List::create(
    Named("param") = NumericVector::create(new_omega, new_rate),
    Named("pdiff") = NumericVector::create(new_omega - omega, new_rate - rate),
    Named("llf") = llf,
    Named("total") = total
  );
}
