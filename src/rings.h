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

/* Stores the projection of x onto y into z. */
void
project(const double x[3], const double y[3], double z[3]);

/* Stores the component of x orthogonal to y in z. */
void
orthogonal_project(const double x[3], const double y[3], double z[3]);

/* Store in y the vector x rotated about the x-axis by an angle theta. */
void
rotate_x(const double x[3], const double theta, double y[3]);

void
rotate_z(const double x[3], const double theta, double y[3]);

/* body.c */

typedef struct {
  double m; /* Mass, in units where G*Mcentral = 1 */
  double a; /* Semi-major axis. */
  double L[3]; /* Of magnitude sqrt(1-e^2), in direction of Lhat */
  double A[3]; /* Of magnitude e, points to periapse. */
} body;

double
mean_motion(const body *b);

void
E_to_rv(const body *b, const double E, double r[3], double v[3]);

/* Fills in the given body from the mass and orbital elements (see
 body structure for units):
 
 * a: semi-major axis.
 * e: eccentricity
 * I: inclination of orbital plane in degrees (I > 90 means retrograde orbit)
 * Omega: Longitude of ascending node in degrees.
 * omega: Argument of periapse in degrees.

*/
void
init_body_from_elements(body *b,
                        const double m, const double a, const double e, const double I, 
                        const double Omega, const double omega);

/* Sets the orbital elements given a body, b. */
void
elements_from_body(const body *b,
                   double *e, double *I, double *Omega, double *omega);

#endif /* __RINGS_H__ */
