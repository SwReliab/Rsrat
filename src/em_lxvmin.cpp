#include <Rcpp.h>
#include <cmath>

#include "gumbel.h"

using namespace Rcpp;

//' @rdname em
// [[Rcpp::export]]

List em_lxvmin_estep(NumericVector params, List data) {
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

  double nn = 0.0;
  double llf = 0.0;
  double t = 0;
  double prev_Fi = 0;
  for (int i=0; i<dsize; i++) {
    t += time[i];
    double Fi = Revd::pgumbel_min(log(t), loc, scale, true, false);
    if (num[i] != 0) {
      nn += num[i];
      llf += num[i] * log(Fi - prev_Fi) - lgamma(num[i]+1);
    }
    if (type[i] == 1) {
      nn += 1;
      llf += Revd::dgumbel_min(log(t), loc, scale, true) - log(t);
    }
    prev_Fi = Fi;
  }
  llf += nn * log(omega) - omega * prev_Fi;
  const double w1 = omega * (1 - prev_Fi);
  const double total = w1;

  return List::create(
    Named("llf") = llf,
    Named("omega") = nn + w1,
    Named("w1") = w1,
    Named("total") = nn + w1
  );
}

//' @rdname em
// [[Rcpp::export]]

double em_lxvmin_pllf(NumericVector params, List data, double w1) {
  const double loc = params[0];
  const double scale = params[1];
  const int dsize = data["len"];
  NumericVector time = as<NumericVector>(data["time"]);
  NumericVector num = as<NumericVector>(data["fault"]);
  IntegerVector type = as<IntegerVector>(data["type"]);

  double llf = 0;
  double prev_Fi = 0;
  double t = 0;
  for (int i=0; i<dsize; i++) {
    t += time[i];
    double Fi = Revd::pgumbel_min(log(t), loc, scale, true, false);
    if (num[i] != 0) {
      llf += num[i] * log(Fi - prev_Fi);
    }
    if (type[i] == 1) {
      llf += Revd::dgumbel_min(log(t), loc, scale, true);
    }
    prev_Fi = Fi;
  }
  llf += w1 * log(1 - prev_Fi);
  return llf;
}
