#include"rings.h"

#include<math.h>

double
dot(const double x[3], const double y[3]) {
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

double
norm(const double x[3]) {
  return sqrt(dot(x,x));
}

void
cross(const double x[3], const double y[3], double z[3]) {
  z[0] = x[1]*y[2] - y[1]*x[2];
  z[1] = x[2]*y[0] - y[2]*x[0];
  z[2] = x[0]*y[1] - y[0]*x[1];
}

void
vscale(const double s, const double x[3], double y[3]) {
  y[0] = x[0]*s;
  y[1] = x[1]*s;
  y[2] = x[2]*s;
}

void
vadd(const double x[3], const double y[3], double z[3]) {
  z[0] = x[0]+y[0];
  z[1] = x[1]+y[1];
  z[2] = x[2]+y[2];
}

void
vsub(const double x[3], const double y[3], double z[3]) {
  z[0] = x[0]-y[0];
  z[1] = x[1]-y[1];
  z[2] = x[2]-y[2];
}

void
unitize(const double x[3], double y[3]) {
  vscale(1.0/norm(x), x, y);
}

void
project(const double x[3], const double y[3], double z[3]) {
  double xdyhat;
  double yhat[3];
  int i;

  unitize(y, yhat);
  xdyhat = dot(x, yhat);

  for (i = 0; i < 3; i++) {
    z[i] = xdyhat*yhat[i];
  }
}

void
orthogonal_project(const double x[3], const double y[3], double z[3]) {
  double px[3];
  int i;

  project(x,y,px);

  for(i = 0; i < 3; i++) {
    z[i] = x[i] - px[i];
  }
}
