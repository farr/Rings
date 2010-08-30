#ifndef __RINGS_H__
#define __RINGS_H__

/* Vector Utilities (vectors.c) */
double
dot(const double x[3], const double y[3]);

double
norm(const double x[3]);

void
cross(const double x[3], const double y[3], double z[3]);

void
vscale(const double s, const double x[3], double y[3]);

void 
vadd(const double x[3], const double y[3], double z[3]);

void
vsub(const double x[3], const double y[3], double z[3]);

void
unitize(const double x[3], double y[3]);

/* body.c */

typedef struct {
  double m; /* Mass, in units where G*Mcentral = 1 */
  double L[3]; /* Of magnitude sqrt(1-e^2), in direction of Lhat */
  double A[3]; /* Of magnitude e, points to periapse. */
} body;

void
EToRv(const body *b, const double E, double r[3], double v[3]);

#endif /* __RINGS_H__ */
