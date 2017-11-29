#include <Rcpp.h>
#include <cmath>

#include "gumbel.h"

using namespace Rcpp;

//' @rdname em
// [[Rcpp::export]]

List em_lxvmax_emstep(NumericVector params, List data) {
  const double omega = params[0];
  const double loc = params[1];
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
  double y = exp(-(log(t)-loc)/scale);
  double en1 = 0.0;
  double en2 = 0.0;
  double en3 = 0.0;
  double llf = 0.0;

  double g00 = 1.0;
  double g01 = exp(-loc/scale);
  double g02 = 1.0;
  double g10 = 1.0 - exp(-y);
  double g11 = exp(-loc/scale)*(1.0 - (1.0 + y)*exp(-y));
  double g12 = 1.0 - exp(-y) * (1.0 + y * log(y));
  if (x != 0.0) {
    double tmp1 = exp(-y);
    double tmp2 = exp(-loc/scale) - g11;
    double tmp3 = exp(-y) * (1.0 + y * log(y));
    en1 += x;
    en2 += x * tmp2 / tmp1;
    en3 += x * tmp3 / tmp1;
    llf += x * log(tmp1) - lgamma(x+1);
  }
  if (type[0] == 1) {
    en1 += 1.0;
    en2 += exp(-log(t)/scale);
    en3 += (log(t)-loc)/scale * (1.0-y);
    llf += Revd::dgumbel(log(t), loc, scale, true) - log(t); // log(lxvmax_pdf(t, scale, loc));
  }
  for (int j=1; j<dsize; j++) {
    x = num[j];
    if (time[j] != 0.0) {
      t += time[j];
      y = exp(-(log(t)-loc)/scale);
      g00 = g10;
      g01 = g11;
      g02 = g12;
      g10 = 1.0 - exp(-y);
      g11 = exp(-loc/scale)*(1.0 - (1.0 + y)*exp(-y));
      g12 = 1.0 - exp(-y) * (1.0 + y * log(y));
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
      en2 += exp(-log(t)/scale);
      en3 += (log(t)-loc)/scale * (1.0-y);
      llf += Revd::dgumbel(log(t), loc, scale, true) - log(t); //log(lxvmax_pdf(t, scale, loc));
    }
  }
  llf += log(omega) * en1;  // en1 is total number of faults
  en1 += omega * g10;  // g10 is the last time
  en2 += omega * g11;  // g11 is the last time
  en3 += omega * g12;  // g12 is the last time
  llf += - omega * (1.0 - g10);

  // M-step
  double total = en1;
  double new_omega =  en1;
  double new_loc = -scale * log(en2/en1);
  //	scale = scale + Math.log(en3) - Math.log(en1);
  double new_scale = scale * en3/en1;

  return List::create(
    Named("param") = NumericVector::create(new_omega, new_loc, new_scale),
    Named("pdiff") = NumericVector::create(new_omega - omega, new_loc - loc, new_scale - scale),
    Named("llf") = llf,
    Named("total") = total
  );
}
