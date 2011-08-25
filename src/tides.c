/* Tidal equations taken from Barker and Ogilvie, On the Tidal
   Evolution of Hot Jupiters on Inclined Orbits, MNRAS, Vol. 395,
   Issue 4, pp 2268--2287, 2009, arXiv:0902.4563, DOI
   10.1111/j.1365-2966.2009.14694.x  */

#include"rings.h"
#include<math.h>

static double pow_int(const double x, const int n) {
  if (n < 0) {
    return pow_int(1.0/x, -n);    
  } else if (n == 0) {
    return 0;
  } else if (n == 1) {
    return x;
  } else if (n % 2 == 0) {
    double sqrt_xn = pow_int(x, n/2);

    return sqrt_xn*sqrt_xn;
  } else {
    double sqrt_xn = pow_int(x, (n-1)/2);

    return x * sqrt_xn * sqrt_xn;
  }
}

static double f1(const double e2) {
  double e4 = e2*e2;
  double e6 = e4*e2;
  double sqrtDenom = sqrt(1.0 - e2);
  double denom = pow_int(sqrtDenom, 13);

  return (1.0 + 15.0/4.0*e2 + 15.0/8.0*e4 + 5.0/64.0*e6)/denom;
}

static double f2(const double e2) {
  double e4 = e2*e2;
  double denom = pow_int(1.0 - e2, 5);

  return (1.0 + 3.0/2.0*e2 + 1.0/8.0*e4)/denom;
}

static double f3(const double e2) {
  double e4 = e2*e2;
  double denom = pow_int(1.0 - e2, 5);

  return (1.0 + 9.0/2.0*e2 + 5.0/8.0*e4)/denom;
}

static double f4(const double e2) {
  double e4 = e2*e2;
  double e6 = e2*e4;
  double denom = pow_int(sqrt(1.0-e2), 13);
  
  return (1.0 + 15.0/2.0*e2 + 45.0/8.0*e4 + 5.0/16.0*e6)/denom;
}

static double f5(const double e2) {
  double denom = pow_int(1.0 - e2, 5);

  return (3.0 + 0.5*e2)/denom;
}

static double f6(const double e2) {
  double e4 = e2*e2;
  double e6 = e4*e2;
  double e8 = e6*e2;
  double denom = pow_int(1.0-e2, 8);

  return (1.0 + 31.0/2.0*e2 + 255.0/8.0*e4 + 185.0/16.0*e6 + 25.0/64.0*e8)/denom;
}
