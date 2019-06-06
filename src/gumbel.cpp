#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

namespace Revd {

  double dgumbel(double x, double loc, double scale, bool log) {
    if (log == false) {
      double y = exp(-(x-loc)/scale);
      return y * exp(-y) / scale;
    } else { // if (log == true) {
      double y = exp(-(x-loc)/scale);
      return -(x-loc)/scale - y - std::log(scale);
    }
  }

  double dgumbel_min(double x, double loc, double scale, bool log) {
    return dgumbel(-x, loc, scale, log);
  }

  double pgumbel(double q, double loc, double scale, bool lower, bool log) {
    if (log == false && lower == true) {
      double y = exp(-(q-loc)/scale);
      return exp(-y);
    } else if (log == false && lower == false) {
      double y = exp(-(q-loc)/scale);
      return 1 - exp(-y);
    } else if (log == true && lower == true) {
      double y = exp(-(q-loc)/scale);
      return -y;
    } else { // (log == true && lower == false) {
      double y = exp(-(q-loc)/scale);
      return std::log(1 - exp(-y));
    }
  }

  double pgumbel_min(double q, double loc, double scale, bool lower, bool log) {
    return pgumbel(-q, loc, scale, !lower, log);
  }

  double qgumbel(double p, double loc, double scale, bool lower, bool log) {
    if (log == true) {
      p = exp(p);
    }
    if (lower == false) {
      p = 1 - p;
    }
    return -scale * std::log(-std::log(p)) + loc;
  }

  double qgumbel_min(double p, double loc, double scale, bool lower, bool log) {
    if (log == true) {
      p = exp(p);
    }
    if (lower == true) {
      p = 1 - p;
    }
    return scale * std::log(-std::log(p)) - loc;
  }
}

//' @rdname gumbel
//' @export
// [[Rcpp::export]]

NumericVector dgumbel(NumericVector x, double loc = 0, double scale = 1, bool log = false, bool min = false) {
  NumericVector result(x.length());
  if (min == false) {
    for (int i=0; i<x.length(); i++) {
      result[i] = Revd::dgumbel(x[i], loc, scale, log);
    }
  } else {
    for (int i=0; i<x.length(); i++) {
      result[i] = Revd::dgumbel_min(x[i], loc, scale, log);
    }
  }
  return result;
}

//' @rdname gumbel
//' @export
// [[Rcpp::export]]

NumericVector pgumbel(NumericVector q, double loc = 0, double scale = 1, bool lower = true, bool log = false, bool min = false) {
  NumericVector result(q.length());
  if (min == false) {
    for (int i=0; i<q.length(); i++) {
      result[i] = Revd::pgumbel(q[i], loc, scale, lower, log);
    }
  } else {
    for (int i=0; i<q.length(); i++) {
      result[i] = Revd::pgumbel_min(q[i], loc, scale, lower, log);
    }
  }
  return result;
}

//' @rdname gumbel
//' @export
// [[Rcpp::export]]

NumericVector qgumbel(NumericVector p, double loc = 0, double scale = 1, bool lower = true, bool log = false, bool min = false) {
  NumericVector result(p.length());
  if (min == false) {
    for (int i=0; i<p.length(); i++) {
      result[i] = Revd::qgumbel(p[i], loc, scale, lower, log);
    }
  } else {
    for (int i=0; i<p.length(); i++) {
      result[i] = Revd::qgumbel_min(p[i], loc, scale, lower, log);
    }
  }
  return result;
}

//' @rdname gumbel
//' @export
// [[Rcpp::export]]

NumericVector rgumbel(int n, double loc = 0, double scale = 1, bool min = false) {
  NumericVector result(n);
  if (min == false) {
    for (int i=0; i<n; i++) {
      double u = R::runif(0,1);
      result[i] = Revd::qgumbel(u, loc, scale, true, false);
    }
  } else {
    for (int i=0; i<n; i++) {
      double u = R::runif(0,1);
      result[i] = Revd::qgumbel_min(u, loc, scale, true, false);
    }
  }
  return result;
}
