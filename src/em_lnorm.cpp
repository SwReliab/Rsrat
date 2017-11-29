#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

//' @rdname em
// [[Rcpp::export]]

List em_lnorm_emstep(NumericVector params, List data) {
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

  double t = time[0];
  double x = num[0];
  double y = (log(t)-mu)/sig;
  double en1 = 0.0;
  double en2 = 0.0;
  double en3 = 0.0;
  double llf = 0.0;
  double tmp = R::dnorm(y, 0, 1, false);
  double g00 = 1.0;
  double g01 = mu;
  double g02 = sig*sig + mu*mu;
  double g10 = R::pnorm(y, 0, 1, false, false);
  double g11 = sig*tmp + mu*g10;
  double g12 = (sig*log(t) + mu*sig)*tmp + (sig*sig + mu*mu)*g10;

  if (x != 0.0) {
    double tmp1 = g00 - g10;
    double tmp2 = g01 - g11;
    double tmp3 = g02 - g12;
    en1 += x;
    en2 += x * tmp2 / tmp1;
    en3 += x * tmp3 / tmp1;
    llf += x * log(tmp1) - lgamma(x+1);
  }
  if (type[0] == 1) {
    en1 += 1.0;
    en2 += log(t);
    en3 += log(t) * log(t);
    llf += R::dnorm(log(t), mu, sig, true) - log(t); // log(lnorm_pdf(t, sig, mu));
  }
  for (int j=1; j<dsize; j++) {
    x = num[j];
    if (time[j] != 0.0) {
      t += time[j];
      y = (log(t)-mu)/sig;
      tmp = R::dnorm(y, 0, 1, false);
      g00 = g10;
      g01 = g11;
      g02 = g12;
      g10 = R::pnorm(y, 0, 1, false, false); // q_normal(y);
      g11 = sig*tmp + mu*g10;
      g12 = (sig*log(t) + mu*sig)*tmp + (sig*sig + mu*mu)*g10;
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
      en2 += log(t);
      en3 += log(t) * log(t);
      llf += R::dnorm(log(t), mu, sig, true) - log(t); // log(lnorm_pdf(t, sig, mu));
    }
  }
  llf += log(omega) * en1;  // en1 is total number of faults
  en1 += omega * g10;  // g10 is the last time
  en2 += omega * g11;  // g11 is the last time
  en3 += omega * g12;  // g12 is the last time
  llf += - omega * (1.0 - g10);

  double total = en1;
  double new_omega =  en1;
  double new_mu = en2 / en1;
  double new_sig = sqrt(en3 / en1 - new_mu*new_mu);

  return List::create(
    Named("param") = NumericVector::create(new_omega, new_mu, new_sig),
    Named("pdiff") = NumericVector::create(new_omega - omega, new_mu - mu, new_sig - sig),
    Named("llf") = llf,
    Named("total") = total
  );
}
