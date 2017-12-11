#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

//
// Marshall-Olkin-type
//

namespace tlogis {
  double func_Fi(double t, double p, double b) {
    return p * (1 - exp(-b*t)) / (p + (1-p) * exp(-b*t));
  }

  double func_h1i(double t, double p, double b) {
    return 2 * (1-p) * (1 - exp(-b*t)) / (p + (1-p) * exp(-b*t));
  }

  double func_h2i(double t, double p, double b) {
    return 2 * (1-p) * (1 - (1+b*t) * exp(-b*t)) / (p + (1-p) * exp(-b*t)) / b;
  }

  double func_H1i(double t, double p, double b) {
    return p * (1 - exp(-b*t)) / (p + (1-p) * exp(-b*t)) / (p + (1-p) * exp(-b*t));
  }

  double func_H2i(double t, double p, double b) {
    return p * (1 - (1+b*t) * exp(-b*t)) / (p + (1-p) * exp(-b*t)) / (p + (1-p) * exp(-b*t)) / b;
  }
}

//' @rdname em
//' @details
//' \code{em_tlogis_emstep_mo} is an emstep based on Marshall-Olkin-type (maximum) with Exp
// [[Rcpp::export]]

List em_tlogis_emstep_mo(NumericVector params, List data) {
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

  const double p = 1/(1 + exp(loc/scale));
  const double b = 1/scale;

  double en1 = 0.0;
  double en2 = 0.0;
  double en3 = 0.0;
  double llf = 0.0;
  double t = 0;
  double prev_Fi = 0;
  double prev_H1i = 0;
  double prev_H2i = 0;
  for (int i=0; i<dsize; i++) {
    t += time[i];
    double Fi = tlogis::func_Fi(t, p, b);
    double H1i = tlogis::func_H1i(t, p, b);
    double H2i = tlogis::func_H2i(t, p, b);

    double x = num[i];
    if (x != 0) {
      double tmp1 = Fi - prev_Fi;
      double tmp2 = H1i - prev_H1i;
      double tmp3 = H2i - prev_H2i;
      en1 += x;
      en2 += x * tmp2 / tmp1;
      en3 += x * tmp3 / tmp1;
      llf += x * log(tmp1) - lgamma(x+1);
    }

    if (type[i] == 1) {
      en1 += 1;
      en2 += 1 + tlogis::func_h1i(t, p, b);
      en3 += 1 + tlogis::func_h2i(t, p, b);
      llf += R::dlogis(t, loc, scale, true) - R::plogis(0, loc, scale, false, true);
    }
    prev_Fi = Fi;
    prev_H1i = H1i;
    prev_H2i = H2i;
  }
  llf += en1 * log(omega) - omega * prev_Fi;
  en1 += omega * (1 - prev_Fi);
  en2 += omega * (1/p - prev_H1i);
  en3 += omega * (1/p/b - prev_H2i);

  const double new_p = en1 / en2;
  const double new_b = en2 / en3;

  const double total = en1;
  const double new_omega = en1;
  const double new_loc = log(1/new_p - 1) / new_b;
  const double new_scale = 1/new_b;

  return List::create(
    Named("param") = NumericVector::create(new_omega, new_loc, new_scale),
    Named("pdiff") = NumericVector::create(new_omega - omega, new_loc - loc, new_scale - scale),
    Named("llf") = llf,
    Named("total") = total
  );
}

//' @rdname em
// [[Rcpp::export]]

List em_tlogis_estep(NumericVector params, List data) {
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

  const double F0 = R::plogis(0, loc, scale, true, false);
  const double barF0 = R::plogis(0, loc, scale, false, false);
  const double log_barF0 = R::plogis(0, loc, scale, false, true);

  double nn = 0.0;
  double llf = 0.0;
  double t = 0;
  double prev_Fi = F0;
  for (int i=0; i<dsize; i++) {
    t += time[i];
    double Fi = R::plogis(t, loc, scale, true, false);
    double x = num[i];
    if (x != 0) {
      nn += x;
      llf += x * (log(Fi - prev_Fi) - log_barF0) - lgamma(x+1);
    }
    if (type[i] == 1) {
      nn += 1;
      llf += R::dlogis(t, loc, scale, true) - log_barF0;
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

double em_tlogis_pllf(NumericVector params, List data, double w0, double w1) {
  const double loc = params[0];
  const double scale = params[1];
  const int dsize = data["len"];
  NumericVector time = as<NumericVector>(data["time"]);
  NumericVector num = as<NumericVector>(data["fault"]);
  IntegerVector type = as<IntegerVector>(data["type"]);

  double llf = w0 * R::plogis(0, loc, scale, true, true);
  double prev = R::plogis(0, loc, scale, true, false);
  double t = 0.0;
  for (int i=0; i<dsize; i++) {
    t += time[i];
    double cur = R::plogis(t, loc, scale, true, false);
    if (num[i] != 0) {
      llf += num[i] * log(cur - prev);
    }
    if (type[i] == 1) {
      llf += R::dlogis(t, loc, scale, true);
    }
    prev = cur;
  }
  llf += w1 * R::plogis(t, loc, scale, false, true);
  return llf;
}
