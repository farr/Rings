#ifndef __RINGS_H__
#define __RINGS_H__

#include<stdlib.h>
#include<gsl/gsl_integration.h>

/* Vector Utilities (vectors.c) */
double
dot(const double x[3], const double y[3]);

double
norm(const double x[3]);

double
distance_squared(const double x[3], const double y[3]);

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

/* Given a vector and its time derivative, returns the change in
   magnitude of the vector. */
double
vdot_to_vmagdot(const double v[3], const double vdot[3]);

/* Returns the acceleration on a body of m = 1 at r1 due to a body of
   m = 1 at r2.  eps is the softening parameter. */
void
softened_specific_acceleration(const double eps, 
                               const double r1[3], const double r2[3],
                               double a[3]);

/* body.c */

typedef struct {
  double m; /* Mass, in units where G*Mcentral = 1 */
  double a; /* Semi-major axis. */
  double L[3]; /* Of magnitude sqrt(1-e^2), in direction of Lhat */
  double A[3]; /* Of magnitude e, points to periapse. */
} body;

#define BODY_VECTOR_SIZE 8
#define BODY_M_INDEX 0
#define BODY_a_INDEX 1
#define BODY_L_INDEX 2
#define BODY_A_INDEX 5

void
body_to_vector(const body *b, double *v);

void
vector_to_body(const double *v, body *b);

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

void
body_instantaneous_rhs(const double eps,
                       const body *b1, const double E1,
                       const body *b2, const double E2,
                       double dbdt[BODY_VECTOR_SIZE]);


/* raw_average.c */

void
raw_average_rhs(const double eps, const body *b1, const body *b2,
                gsl_integration_workspace *ws1, const size_t ws1_limit,
                gsl_integration_workspace *ws2, const size_t ws2_limit,
                double rhs[BODY_VECTOR_SIZE]);

#endif /* __RINGS_H__ */
