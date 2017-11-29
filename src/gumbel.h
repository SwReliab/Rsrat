#ifndef _GUMBEL_H_
#define _GUMBEL_H_

namespace Revd {

  double dgumbel(double x, double loc, double scale, bool log);
  double dgumbel_min(double x, double loc, double scale, bool log);
  double pgumbel(double q, double loc, double scale, bool lower, bool log);
  double pgumbel_min(double q, double loc, double scale, bool lower, bool log);
  double qgumbel(double p, double loc, double scale, bool lower, bool log);
  double qgumbel_min(double p, double loc, double scale, bool lower, bool log);

}

#endif
