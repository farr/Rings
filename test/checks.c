#include"checks.h"
#include<math.h>

int
check_close(const double epsabs, const double epsrel, const double x, const double y) {
  double dx = fabs(x-y);
  double ave = 0.5*(fabs(x)+fabs(y));

  return dx <= epsabs + epsrel*ave;
}

int
check_vector_close(const double epsabs, const double epsrel, const double x[3], const double y[3]) {
  int i;
  for (i = 0; i < 3; i++) {
    if (!check_close(epsabs, epsrel, x[i], y[i])) return 0;
  }

  return 1;
}
