// inverse function of log(a) - psi(a)

#include "math.h"
#include "numlib.h"

#define PI (3.14159265358979324)
#define GAM (0.577215664901532861)

#define EPS 1.0e-8
#define MAXCNT 200

double inv_log_psi(double v) {
  double init;
  if (v > 1.0) {
    // approximation form 1/x + (x-1)
    init = (1.0 + v - sqrt(v*v + 2.0*v - 3.0))/2.0;
  } else {
    // approximation form GAM / (1 - (PI^2 - 6) (x-1) / (6 GAM))
    init = (6.0 * GAM * GAM - 6.0 * v - 6.0 * GAM * v + PI * PI * v) / ((PI * PI - 6.0)*v);
  }
  int cnt;
  double a, b, c;
  a = init/2.0;
  b = init;
  cnt = 0;
  while (log(b) - psi(b) > v
	 && cnt++ < MAXCNT) {
    a = b;
    b *= 2.0;
  }
  while (log(a) - psi(a) <= v
	 && cnt++ < MAXCNT) {
    b = a;
    a = a/2.0;
  }
  cnt = 0;
  c = (a+b)/2.0;
  while (fabs(a - b)/a > EPS
	 && cnt++ < MAXCNT) {
    c = (a+b)/2.0;
    if (log(c) - psi(c) < v) {
      b = c;
    } else {
      a = c;
    }
  }
  return c;
}
