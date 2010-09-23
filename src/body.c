#include"rings.h"
#include<assert.h>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>

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

/* body * */
/* alloc_body_from_elements(const double m, const double a, const double e, const double I, */
/*                          const double Omega, const double omega) { */
/*   double Lmag = sqrt(1.0-e*e); */
/*   body *res = malloc(sizeof(body)); */

/*   assert(res != 0); */

/*   res->m = m; */
/*   res->a = a; */

/*   res->L[2] = cos(I*M_PI/180.0); */
/*   res->L[0] = cos(Omega*M_PI/180.0 - M_PI/2); */
/*   res->L[1] = sin(Omega*M_PI/180.0 - M_PI/2); */

/*   unitize(res->L); */
/*   vscale(Lmag, res->L, res->L); */

  
/* } */

/* void */
/* elements_from_body(const body *b, */
/*                    double *e, double *I, double *Omega, double *omega) { */
  
/*   double Lmag = norm(b->L); */
/*   double lang = atan2(b->L[1], b->L[0]); */
/*   double asc_node[3]; */
/*   double my_omega; */
/*   double asc_node_x_A[3]; */

/*   *I = acos(b->L[2]/Lmag)*180.0/M_PI; */
/*   *e = norm(b->A); */
/*   *Omega = fmod(lang+M_PI/2.0, 2.0*M_PI)*180.0/M_PI; */
  
/*   my_omega = *omega; */

/*   asc_node[0] = cos(my_omega); */
/*   asc_node[1] = sin(my_omega); */
/*   asc_node[2] = 0.0; */

/*   *omega = acos(dot(asc_node, b->A)/norm(b->A)); */

/*   cross(asc_node, b->A, asc_node_x_A); */
/*   if (dot(b->L, asc_node_x_A) < 0) { /\* This means large argument of periapse. *\/ */
/*     *omega = 2.0*M_PI - *omega; */
/*   } */

/*   *omega = *omega*180.0/M_PI; */
/* } */
