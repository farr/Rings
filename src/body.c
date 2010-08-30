#include"rings.h"
#include<assert.h>
#include<math.h>

double
meanMotion(const body *b) {
  double a = b->a;
  return sqrt((1.0+b->m)/(a*a*a));
}

void
EtoRv(const body *b, const double E, double x[3], double v[3]) {
  double xhat[3], yhat[3], zhat[3];
  double e = norm(b->A);
  double n = meanMotion(b);
  double cosE = cos(E), sinE = sin(E);
  double eCosE = e*cosE;
  double sqrt1me2 = norm(b->L);
  double denom = 1.0 - eCosE;
  double a = b->a;
  double xc = a*(cosE - e);
  double yc = a*sqrt1me2*sinE;
  double r = a*(1.0 - eCosE);
  double na2 = n*a*a;
  double xdot = -na2*sinE/r;
  double ydot = na2*sqrt1me2*cosE/r;

  int i;

  assert(e > 0.0);
  assert(dot(b->A, b->L)/e/sqrt1me2 < 1e-8);

  unitize(b->L, zhat);
  unitize(b->A, xhat);
  cross(zhat, xhat, yhat);

  for (i = 0; i < 3; i++) {
    x[i] = xc*xhat[i] + yc*yhat[i];
    v[i] = xdot*xhat[i] + ydot*yhat[i];
  }
}
