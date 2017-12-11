#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

namespace tnorm {
  double func_barFi(double t, double mu, double sig) {
    return R::pnorm(t, mu, sig, false, false);
  }

  double func_barPhi1i(double t, double mu, double sig) {
    double z = (t - mu) / sig;
    return sig * R::dnorm(z, 0, 1, false) + mu * R::pnorm(z, 0, 1, false, false);
  }

  double func_barPhi2i(double t, double mu, double sig) {
    double z = (t - mu) / sig;
    return (sig*t + sig*mu) * R::dnorm(z, 0, 1, false) + (sig*sig + mu*mu) * R::pnorm(z, 0, 1, false, false);
  }
}

//' @rdname em
//' @details
//' \code{em_tnorm_emstep} has been modified by using EM for truncated distribution. Concretely,
//' when the distribution is truncated at origin, the expected value is given by
//' \deqn{\text{E}[h(X)|D] = \sum_{i=1}^k h(x_i) + \frac{k}{\overline{F}(0)} \int_{-\infty}^0 h(x) f(x) dx.}
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

  const double barF0 = tnorm::func_barFi(0, mu, sig);
  const double barPhi10 = tnorm::func_barPhi1i(0, mu, sig);
  const double barPhi20 = tnorm::func_barPhi2i(0, mu, sig);

  double nn = 0;
  double en1 = 0;
  double en2 = 0;
  double en3 = 0;
  double llf = 0;

  double t = 0;
  double prev_barFi = barF0;
  double prev_barPhi1i = barPhi10;
  double prev_barPhi2i = barPhi20;

  for (int i=0; i<dsize; i++) {
    t += time[i];
    double barFi = tnorm::func_barFi(t, mu, sig);
    double barPhi1i = tnorm::func_barPhi1i(t, mu, sig);
    double barPhi2i = tnorm::func_barPhi2i(t, mu, sig);

    double x = num[i];
    if (x != 0) {
      double tmp1 = prev_barFi - barFi;
      double tmp2 = prev_barPhi1i - barPhi1i;
      double tmp3 = prev_barPhi2i - barPhi2i;
      nn += x;
      en1 += x;
      en2 += x * tmp2 / tmp1;
      en3 += x * tmp3 / tmp1;
      llf += x * log(tmp1) - lgamma(x+1);
    }
    if (type[i] == 1) {
      nn += 1;
      en1 += 1;
      en2 += t;
      en3 += t * t;
      llf += R::dnorm(t, mu, sig, true);
    }
    prev_barFi = barFi;
    prev_barPhi1i = barPhi1i;
    prev_barPhi2i = barPhi2i;
  }
  llf += nn * (log(omega) - log(barF0)) - omega * (barF0 - prev_barFi) / barF0;
  en1 += omega * prev_barFi / barF0;
  en2 += omega * prev_barPhi1i / barF0;
  en3 += omega * prev_barPhi2i / barF0;

  double en1dash = en1 + en1 * (1 - barF0) / barF0;
  en2 += en1 * (mu - barPhi10) / barF0;
  en3 += en1 * (sig*sig + mu*mu - barPhi20) / barF0;

  double new_mu = en2 / en1dash;
  double new_sig = sqrt(en3 / en1dash - new_mu*new_mu);
  double new_omega = en1;
  double total = en1;

  return List::create(
    Named("param") = NumericVector::create(new_omega, new_mu, new_sig),
    Named("pdiff") = NumericVector::create(new_omega - omega, new_mu - mu, new_sig - sig),
    Named("llf") = llf,
    Named("total") = total
  );
}
