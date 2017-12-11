#include <Rcpp.h>
#include <cmath>

#include "gumbel.h"

using namespace Rcpp;

//' @rdname em
// [[Rcpp::export]]

List em_txvmin_estep(NumericVector params, List data) {
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

  const double F0 = Revd::pgumbel_min(0, loc, scale, true, false);
  const double barF0 = Revd::pgumbel_min(0, loc, scale, false, false);
  const double log_barF0 = Revd::pgumbel_min(0, loc, scale, false, true);

  double nn = 0.0;
  double llf = 0.0;
  double t = 0;
  double prev_Fi = F0;
  for (int i=0; i<dsize; i++) {
    t += time[i];
    double Fi = Revd::pgumbel_min(t, loc, scale, true, false);
    double x = num[i];
    if (x != 0) {
      nn += x;
      llf += x * (log(Fi - prev_Fi) - log_barF0) - lgamma(x+1);
    }
    if (type[i] == 1) {
      nn += 1;
      llf += Revd::dgumbel_min(t, loc, scale, true) - log_barF0;
    }
    prev_Fi = Fi;
  }
  llf += nn * log(omega) - omega * (prev_Fi - F0) / barF0;
  const double w1 = omega * (1 - prev_Fi) / barF0;
  const double w0 = (nn + w1) * F0 / barF0;

  return List::create(
    Named("llf") = llf,
    Named("omega") = nn + w1,
    Named("w0") = w0,
    Named("w1") = w1,
    Named("total") = nn + w1
  );
}

//' @rdname em
// [[Rcpp::export]]

double em_txvmin_pllf(NumericVector params, List data, double w0, double w1) {
  const double loc = params[0];
  const double scale = params[1];
  const int dsize = data["len"];
  NumericVector time = as<NumericVector>(data["time"]);
  NumericVector num = as<NumericVector>(data["fault"]);
  IntegerVector type = as<IntegerVector>(data["type"]);

  double llf = w0 * Revd::pgumbel_min(0, loc, scale, true, true);
  double prev = Revd::pgumbel_min(0, loc, scale, true, false);
  double t = 0.0;
  for (int i=0; i<dsize; i++) {
    t += time[i];
    double cur = Revd::pgumbel_min(t, loc, scale, true, false);
    if (num[i] != 0) {
      llf += num[i] * log(cur - prev);
    }
    if (type[i] == 1) {
      llf += Revd::dgumbel_min(t, loc, scale, true);
    }
    prev = cur;
  }
  llf += w1 * Revd::pgumbel_min(t, loc, scale, false, true);
  return llf;
}
