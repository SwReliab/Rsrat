#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

//' @rdname em
// [[Rcpp::export]]

List em_tnorm_emstep(NumericVector params, List data) {
  const double omega = params[0];
  const double mu = params[1];
  const double sig = params[2];
  const int dsize = data["len"];
  NumericVector time = as<NumericVector>(data["time"]);
  NumericVector num = as<NumericVector>(data["fault"]);
  IntegerVector type = as<IntegerVector>(data["type"]);

  if (dsize != time.length() || dsize != num.length() || dsize != type.length()) {
    stop("Invalid data.");
  }

  const double omega_dash = omega / R::pnorm(0, mu, sig, false, false);

  double en1 = 0.0;
  double en2 = 0.0;
  double en3 = 0.0;
  double llf = 0.0;
  double y = -mu/sig;
  double tmp = R::dnorm(y, 0, 1, false);
  double g00 = R::pnorm(y, 0, 1, false, false);
  double g01 = sig*tmp + mu*g00;
  double g02 = (mu*sig)*tmp + (sig*sig + mu*mu)*g00;

  double t = time[0];
  double x = num[0];
  y = (t - mu) / sig;
  tmp = R::dnorm(y, 0, 1, false);
  double g10 = g00;
  double g11 = g01;
  double g12 = g02;
  double g20 = R::pnorm(y, 0, 1, false, false);
  double g21 = sig*tmp + mu*g20;
  double g22 = (sig*t + mu*sig)*tmp + (sig*sig + mu*mu)*g20;

  if (x != 0.0) {
    double tmp1 = g10 - g20;
    double tmp2 = g11 - g21;
    double tmp3 = g12 - g22;
    en1 += x;
    en2 += x * tmp2 / tmp1;
    en3 += x * tmp3 / tmp1;
    llf += x * (log(tmp1) - log(g00)) - lgamma(x+1.0);
  }
  if (type[0] == 1) {
    en1 += 1.0;
    en2 += t;
    en3 += t * t;
    llf += R::dnorm(t, mu, sig, true) - R::pnorm(0, mu, sig, false, true); // log(tnorm_pdf(t, sig, mu));
  }
  for (int j=1; j<dsize; j++) {
    x = num[j];
    if (time[j] != 0.0) {
      t += time[j];
      y = (t - mu) / sig;
      tmp = R::dnorm(y, 0, 1, false);
      g10 = g20;
      g11 = g21;
      g12 = g22;
      g20 = R::pnorm(y, 0, 1, false, false);
      g21 = sig*tmp + mu*g20;
      g22 = (sig*t + mu*sig)*tmp + (sig*sig + mu*mu)*g20;
    }
    if (x != 0.0) {
      double tmp1 = g10 - g20;
      double tmp2 = g11 - g21;
      double tmp3 = g12 - g22;
      en1 += x;
      en2 += x * tmp2 / tmp1;
      en3 += x * tmp3 / tmp1;
      llf += x * (log(tmp1) - log(g00))- lgamma(x+1.0);
    }
    if (type[j] == 1) {
      en1 += 1.0;
      en2 += t;
      en3 += t * t;
      llf += R::dnorm(t, mu, sig, true) - R::pnorm(0, mu, sig, false, true); // log(tnorm_pdf(t, sig, mu));
    }
  }
  llf += (log(omega_dash) + log(g00)) * en1;  // en1 is total number of faults
  double total = en1 + omega_dash * g20;
  en1 += omega_dash * (1.0 - g00 + g20);  // g00 is the first, g10 is the last
  en2 += omega_dash * (mu - g01 + g21);  // g01 is the first, g11 is the last
  en3 += omega_dash * (sig*sig + mu*mu - g02 + g22);  // g02 is the first, g22 is the last
  llf += - omega_dash * (g00 - g20);

  double new_mu = en2 / en1;
  double new_sig = sqrt(en3 / en1 - new_mu*new_mu);
  double new_omega =  en1 * R::pnorm(0, new_mu, new_sig, false, false); //q_normal(-mu/sig);

  return List::create(
    Named("param") = NumericVector::create(new_omega, new_mu, new_sig),
    Named("pdiff") = NumericVector::create(new_omega - omega, new_mu - mu, new_sig - sig),
    Named("llf") = llf,
    Named("total") = total
  );
}
