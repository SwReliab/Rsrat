#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

//' @rdname em
// [[Rcpp::export]]

List em_llogis_emstep(NumericVector params, List data) {
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
  double y = exp((log(t)-loc)/scale);
  double en1 = 0.0;
  double en2 = 0.0;
  double en3 = 0.0;
  double llf = 0.0;
  double g00 = 1.0;
  double g01 = 0.5;
  double g02 = 1.0;
  double g10 = 1.0/(1.0+y);
  double g11 = 1.0/(2.0*(1.0+y)*(1.0+y));
  double g12 = (1.0+(1.0+log(y))*y)/((1.0+y)*(1.0+y));

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
    en2 += 1.0/(1.0+y);
    en3 += (y-1.0)*log(y)/(1.0+y);
    llf += R::dlogis(log(t), loc, scale, true) - log(t); // log(llogist_pdf(t, scale, loc));
  }
  for (int j=1; j<dsize; j++) {
    x = num[j];
    if (time[j] != 0.0) {
      t += time[j];
      y = exp((log(t)-loc)/scale);
      g00 = g10;
      g01 = g11;
      g02 = g12;
      g10 = 1.0/(1.0+y);
      g11 = 1.0/(2.0*(1.0+y)*(1.0+y));
      g12 = (1.0+(1.0+log(y))*y)/((1.0+y)*(1.0+y));
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
      en2 += 1.0/(1.0+y);
      en3 += (y-1.0)*log(y)/(1.0+y);
      llf += R::dlogis(log(t), loc, scale, true) - log(t); // log(llogist_pdf(t, scale, loc));
    }
  }
  llf += log(omega) * en1;  // en1 is total number of faults
  en1 += omega * g10;  // g00 is the last time
  en2 += omega * g11;  // g01 is the last time
  en3 += omega * g12;  // g02 is the last time
  llf += - omega * (1.0 - g10);

  // M-step
  double total = en1;
  double new_omega =  en1;
  double new_scale = scale*en3/en1;
  double new_loc = loc + new_scale*(log(en1/2.0) - log(en2));

  return List::create(
    Named("param") = NumericVector::create(new_omega, new_loc, new_scale),
    Named("pdiff") = NumericVector::create(new_omega - omega, new_loc - loc, new_scale - scale),
    Named("llf") = llf,
    Named("total") = total,
    Named("residual") = omega * g10
  );
}

//' @rdname em
// [[Rcpp::export]]

List em_llogis_estep(NumericVector params, List data) {
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
    double Fi = R::plogis(log(t), loc, scale, true, false);
    if (num[i] != 0) {
      nn += num[i];
      llf += num[i] * log(Fi - prev_Fi) - lgamma(num[i]+1);
    }
    if (type[i] == 1) {
      nn += 1;
      llf += R::dlogis(log(t), loc, scale, true) - log(t);
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

double em_llogis_pllf(NumericVector params, List data, double w1) {
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
    double Fi = R::plogis(log(t), loc, scale, true, false);
    if (num[i] != 0) {
      llf += num[i] * log(Fi - prev_Fi);
    }
    if (type[i] == 1) {
      llf += R::dlogis(log(t), loc, scale, true);
    }
    prev_Fi = Fi;
  }
  llf += w1 * log(1 - prev_Fi);
  return llf;
}
