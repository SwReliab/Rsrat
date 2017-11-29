#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

//' @rdname em
// [[Rcpp::export]]

List em_tlogis_emstep(NumericVector params, List data) {
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

  const double omega_dash = omega / R::plogis(0, loc, scale, false, false);

  double y0 = exp(-loc/scale);
  double en1 = 0.0;
  double en2 = 0.0;
  double en3 = 0.0;
  double llf = 0.0;
  double g00 = 1.0/(1.0+y0);
  double g01 = 1.0/(2.0*(1.0+y0)*(1.0+y0));
  double g02 = (1.0+(1.0+log(y0))*y0)/((1.0+y0)*(1.0+y0));

  double t = time[0];
  double x = num[0];
  double y = exp((t-loc)/scale);
  double g10 = g00;
  double g11 = g01;
  double g12 = g02;
  double g20 = 1.0/(1.0+y);
  double g21 = 1.0/(2.0*(1.0+y)*(1.0+y));
  double g22 = (1.0+(1.0+log(y))*y)/((1.0+y)*(1.0+y));
  if (x != 0.0) {
    double tmp1 = g10 - g20;
    double tmp2 = g11 - g21;
    double tmp3 = g12 - g22;
    en1 += x;
    en2 += x * tmp2 / tmp1;
    en3 += x * tmp3 / tmp1;
    llf += x * (log(tmp1) - log(g00)) - lgamma(x+1);
  }
  if (type[0] == 1) {
    en1 += 1.0;
    en2 += 1.0/(1.0+y);
    en3 += (y-1.0)*log(y)/(1.0+y);
    llf += R::dlogis(t, loc, scale, true) - R::plogis(0, loc, scale, false, true); // log(tlogist_pdf(t, scale, loc));
  }
  for (int j=1; j<dsize; j++) {
    x = num[j];
    if (time[j] != 0.0) {
      t += time[j];
      y = exp((t-loc)/scale);
      g10 = g20;
      g11 = g21;
      g12 = g22;
      g20 = 1.0/(1.0+y);
      g21 = 1.0/(2.0*(1.0+y)*(1.0+y));
      g22 = (1.0+(1.0+log(y))*y)/((1.0+y)*(1.0+y));
    }
    if (x != 0.0) {
      double tmp1 = g10 - g20;
      double tmp2 = g11 - g21;
      double tmp3 = g12 - g22;
      en1 += x;
      en2 += x * tmp2 / tmp1;
      en3 += x * tmp3 / tmp1;
      llf += x * (log(tmp1) - log(g00))- lgamma(x+1);
    }
    if (type[j] == 1) {
      en1 += 1.0;
      en2 += 1.0/(1.0+y);
      en3 += (y-1.0)*log(y)/(1.0+y);
      llf += R::dlogis(t, loc, scale, true) - R::plogis(0, loc, scale, false, true); // log(tlogist_pdf(t, scale, loc));
    }
  }
  llf += (log(omega_dash) + log(g00)) * en1;  // en1 is total number of faults
  double total = en1 + omega_dash * g20;
  en1 += omega_dash * (1.0 - g00 + g20);  // g00 is the first, g10 is the last
  en2 += omega_dash * (0.5 - g01 + g21);  // g01 is the first, g11 is the last
  en3 += omega_dash * (1.0 - g02 + g22);  // g02 is the first, g12 is the last
  llf += - omega_dash * (g00 - g20);

  double new_scale = scale*en3/en1;
  double new_loc = loc + new_scale*(log(en1/2.0) - log(en2));
  y0 = exp(-new_loc/new_scale);
  double new_omega = en1 / (1.0 + y0);

  return List::create(
    Named("param") = NumericVector::create(new_omega, new_loc, new_scale),
    Named("pdiff") = NumericVector::create(new_omega - omega, new_loc - loc, new_scale - scale),
    Named("llf") = llf,
    Named("total") = total
  );
}
