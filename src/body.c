#include"rings.h"
#include<assert.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

double
mean_motion(const body *b) {
  double a = b->a;
  double m = b->m;
  double mm = sqrt((1.0+m)/(a*a*a));
  return mm;
}

void
E_to_rv(const body *b, const double E, double x[3], double v[3]) {
  double xhat[3], yhat[3], zhat[3];
  double e = norm(b->A);
  double n = mean_motion(b);
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
  assert(fabs(sqrt1me2 - sqrt(1.0-e*e)) < 1e-8);

  unitize(b->L, zhat);
  unitize(b->A, xhat);
  cross(zhat, xhat, yhat);

  for (i = 0; i < 3; i++) {
    x[i] = xc*xhat[i] + yc*yhat[i];
    v[i] = xdot*xhat[i] + ydot*yhat[i];
  }
}

static void
rotate_to_orbit_frame(double v[3], const double I, const double Omega, const double omega) {
  double temp[3];

  memcpy(temp, v, 3*sizeof(double));

  rotate_z(temp, omega*M_PI/180.0, v); /* Periapse to omega in x-y plane. */
  rotate_x(v, I*M_PI/180.0, temp); /* Inclined orbit; line of nodes still x-axis. */
  rotate_z(temp, Omega*M_PI/180.0, v); /* Rotate line of nodes to correct location. */
  /* Correct vector stored in v. */
}

body *
alloc_body_from_elements(const double m, const double a, const double e, const double I,
                         const double Omega, const double omega) {
  body *res = malloc(sizeof(body));

  assert(res != 0);

  res->m = m;
  res->a = a;
  
  res->L[0] = 0.0;
  res->L[1] = 0.0;
  res->L[2] = sqrt(1.0 - e*e);

  res->A[0] = e;
  res->A[1] = 0.0;
  res->A[2] = 0.0;

  rotate_to_orbit_frame(res->L, I, Omega, omega);
  rotate_to_orbit_frame(res->A, I, Omega, omega);

  return res;
}

void
elements_from_body(const body *b,
                   double *e, double *I, double *Omega, double *omega) {
  
  double Lmag = norm(b->L);
  double lang = atan2(b->L[1], b->L[0]); /* Between -Pi and Pi */
  double asc_node[3];
  double my_omega;
  double asc_node_x_A[3];

  *I = acos(b->L[2]/Lmag)*180.0/M_PI;
  *e = norm(b->A);
  *Omega = lang + M_PI/2.0;

  if (*Omega < 0.0) {
    *Omega = 2.0*M_PI - *Omega;
  }

  my_omega = *Omega;

  *Omega *= 180.0/M_PI;

  asc_node[0] = cos(my_omega);
  asc_node[1] = sin(my_omega);
  asc_node[2] = 0.0;

  *omega = acos(dot(asc_node, b->A)/norm(b->A));

  cross(asc_node, b->A, asc_node_x_A);
  if (dot(b->L, asc_node_x_A) < 0.0) { /* This means large argument of periapse. */
    *omega = 2.0*M_PI - *omega;
  }

  *omega = *omega*180.0/M_PI;
}
