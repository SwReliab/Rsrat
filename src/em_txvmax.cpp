#include <Rcpp.h>
#include <cmath>

#include "gumbel.h"

using namespace Rcpp;

//
// Marshall-Olkin-type
//

namespace txvmax {
  double func_Fi(double t, double nu, double b) {
    const double z = exp(-b*t);
    return (exp(-nu*z) - exp(-nu)) / (1 - exp(-nu));
  }

  double func_H1i(double t, double nu, double b) {
    const double z = exp(-b*t);
    return nu * exp(-nu*z) * (1 - z) / (1 - exp(-nu));
  }

  double func_H2i(double t, double nu, double b) {
    const double z = exp(-b*t);
    return nu * exp(-nu*z) * (1 - (1+b*t)*z) / (1 - exp(-nu)) / b;
  }

  double func_h1i(double t, double nu, double b) {
    const double z = exp(-b*t);
    return nu * (1 - z);
  }

  double func_h2i(double t, double nu, double b) {
    const double z = exp(-b*t);
    return nu * (1 - (1+b*t)*z) / b;
  }

  double phi(double nu) {
    return nu / (1 - exp(-nu));
  }

  double inv_phi(double x) {
    const double epsi = 1.0e-8;
    double left = x - 1;
    double right = x;
    while (fabs(left - right) > epsi) {
      double med = (left + right) / 2;
      if (phi(med) < x) {
        left = med;
      } else {
        right = med;
      }
    }
    return (left + right) / 2;
  }
}

//' @rdname em
//' @details
//' \code{em_txvmax_emstep_mo} is an emstep based on Marshall-Olkin-type (maximum) with Exp
// [[Rcpp::export]]

List em_txvmax_emstep_mo(NumericVector params, List data) {
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

  const double nu = exp(loc/scale);
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
    double Fi = txvmax::func_Fi(t, nu, b);
    double H1i = txvmax::func_H1i(t, nu, b);
    double H2i = txvmax::func_H2i(t, nu, b);

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
      en2 += 1 + txvmax::func_h1i(t, nu, b);
      en3 += 1 + txvmax::func_h2i(t, nu, b);
      llf += Revd::dgumbel(t, loc, scale, true) - Revd::pgumbel(0, loc, scale, false, true);
    }
    prev_Fi = Fi;
    prev_H1i = H1i;
    prev_H2i = H2i;
  }
  llf += en1 * log(omega) - omega * prev_Fi;
  en1 += omega * (1 - prev_Fi);
  en2 += omega * (txvmax::phi(nu) - prev_H1i);
  en3 += omega * (txvmax::phi(nu)/b - prev_H2i);

  const double new_nu = txvmax::inv_phi(en2 / en1);
  const double new_b = en2 / en3;

  const double total = en1;
  const double new_omega = en1;
  const double new_loc = log(new_nu) / new_b;
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

List em_txvmax_estep(NumericVector params, List data) {
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

  const double F0 = Revd::pgumbel(0, loc, scale, true, false);
  const double barF0 = Revd::pgumbel(0, loc, scale, false, false);
  const double log_barF0 = Revd::pgumbel(0, loc, scale, false, true);

  double nn = 0.0;
  double llf = 0.0;
  double t = 0;
  double prev_Fi = F0;
  for (int i=0; i<dsize; i++) {
    t += time[i];
    double Fi = Revd::pgumbel(t, loc, scale, true, false);
    double x = num[i];
    if (x != 0) {
      nn += x;
      llf += x * (log(Fi - prev_Fi) - log_barF0) - lgamma(x+1);
    }
    if (type[i] == 1) {
      nn += 1;
      llf += Revd::dgumbel(t, loc, scale, true) - log_barF0;
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

double em_txvmax_pllf(NumericVector params, List data, double w0, double w1) {
  const double loc = params[0];
  const double scale = params[1];
  const int dsize = data["len"];
  NumericVector time = as<NumericVector>(data["time"]);
  NumericVector num = as<NumericVector>(data["fault"]);
  IntegerVector type = as<IntegerVector>(data["type"]);

  double llf = w0 * Revd::pgumbel(0, loc, scale, true, true);
  double prev = Revd::pgumbel(0, loc, scale, true, false);
  double t = 0.0;
  for (int i=0; i<dsize; i++) {
    t += time[i];
    double cur = Revd::pgumbel(t, loc, scale, true, false);
    if (num[i] != 0) {
      llf += num[i] * log(cur - prev);
    }
    if (type[i] == 1) {
      llf += Revd::dgumbel(t, loc, scale, true);
    }
    prev = cur;
  }
  llf += w1 * Revd::pgumbel(t, loc, scale, false, true);
  return llf;
}
