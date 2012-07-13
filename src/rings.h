#ifndef __RINGS_H__
#define __RINGS_H__

#include<stdlib.h>
#include<gsl/gsl_integration.h>
#include<gsl/gsl_odeiv.h>
#include<stdio.h>
#include<setjmp.h>

/* Vector Utilities (vectors.c) */

/** Dot product. */
double
dot(const double x[3], const double y[3]);

/** Norm. */
double
norm(const double x[3]);

/** Square of distance. */
double
distance_squared(const double x[3], const double y[3]);

/** \f$\mathbf{z} = \mathbf{x} \times \mathbf{y}\f$ */
void
cross(const double x[3], const double y[3], double z[3]);

/** \f$\mathbf{y} = s \mathbf{x}\f$ */
void
vscale(const double s, const double x[3], double y[3]);

/** \f$\mathbf{z} = \mathbf{x} + \mathbf{y}\f$ */
void 
vadd(const double x[3], const double y[3], double z[3]);

/** \f$\mathbf{z} = \mathbf{x} - \mathbf{y}\f$ */
void
vsub(const double x[3], const double y[3], double z[3]);

/** \f$\mathbf{y} = \hat{\mathbf{x}}\f$ */
void
unitize(const double x[3], double y[3]);

/** Stores the projection of x onto y into z: \f$\mathbf{z} =
    \left(\mathbf{x} \cdot \mathbf{y} \right) \hat{\mathbf{y}}\f$. */
void
project(const double x[3], const double y[3], double z[3]);

/** Stores the component of x orthogonal to y in z: \f$\mathbf{z} =
   \mathbf{x} - \left(\mathbf{x}\cdot\hat{\mathbf{y}}\right)
   \hat{\mathbf{y}}\f$ . */
void
orthogonal_project(const double x[3], const double y[3], double z[3]);

/** \f$\mathbf{u} \leftarrow R_{\hat{x}}(\theta) \mathbf{w}\f$ */
void
rotate_x(const double w[3], const double theta, double u[3]);

/** \f$\mathbf{u} \leftarrow R_{\hat{y}}(\theta) \mathbf{w}\f$ */
void
rotate_z(const double w[3], const double theta, double u[3]);

/** Given a vector and its time derivative, returns the change in
   magnitude of the vector. */
double
vdot_to_vmagdot(const double v[3], const double vdot[3]);

/** Returns the acceleration on a body of m = 1 at r1 due to a body of
    m = 1 at r2.  eps is the softening parameter. */
void
softened_specific_acceleration(const double eps, 
                               const double r1[3], const double r2[3],
                               double a[3]);

/* body.c */

/** Bodies. */
typedef struct {
  double m; /** Mass, in units where G*Mcentral = 1 */
  double a; /** Semi-major axis. */
  double Qp; /** Q' from Barker and Ogilvie.  Q' = 3 Q / (2 k), with Q
                 = 2 Pi E / Edot over an orbit and k the love number.
                 A homogeneous solid body, Q' = Q.  We assume
                 (following Barker and Ogilvie) that Q' does not
                 change as the orbit changes (though between
                 integration steps the user can change Q' to alter the
                 tidal behavior of the planet).  This is equivalent to
                 a small time-lag for all tidal components that scales
                 with the orbital period, so Q ~ 1/(tau*omega) remains
                 constant. */
  double I;  /** Body moment of inertia (units are m*a^2). */
  double R; /** Body radius (in units of a). */
  double L[3]; /** Of magnitude sqrt(1-e^2), in direction of Lhat */
  double A[3]; /** Of magnitude e, points to periapse. */
  double spin[3]; /** Body's instantaneous spin vector. */  
} body;

#define BODY_VECTOR_SIZE 14
#define BODY_M_INDEX 0
#define BODY_a_INDEX 1
#define BODY_Qp_INDEX 2
#define BODY_I_INDEX 3
#define BODY_R_INDEX 4
#define BODY_L_INDEX 5
#define BODY_A_INDEX 8
#define BODY_SPIN_INDEX 11

/** Central body information. */
typedef struct {
  double Qp;  /** Q' from Baker and Ogilvie (see above). */
  double I;  /** Moment of inertia (units are m*a^2). */
  double R;  /** Radius (units are a) */
  double spin[3];  /** Instantaneous spin vector. */
} central_body;

#define CENTRAL_BODY_VECTOR_SIZE 6
#define CENTRAL_BODY_Qp_INDEX 0
#define CENTRAL_BODY_I_INDEX 1
#define CENTRAL_BODY_R_INDEX 2
#define CENTRAL_BODY_SPIN_INDEX 3

/** Unpack a body into a vector of length #BODY_VECTOR_SIZE */
void
body_to_vector(const body *b, double *v);

/** Unpack a central body into a vector of length #CENTRAL_BODY_VECTOR_SIZE */
void
central_body_to_vector(const central_body *bc, double *v);

/** Pack a vector of length #BODY_VECTOR_SIZE into b. */
void
vector_to_body(const double *v, body *b);

/** Pack a vector of length #CENTRAL_BODY_VECTOR_SIZE */
void
vector_to_central_body(const double *v, central_body *bc);

/** Mean motion of b. */
double
mean_motion(const body *b);

/** Returns the eccentricity of the given body in an efficient and
    accurate way for both \f$e \to 0\f$ and \f$e \to 1\f$.  Returns
    NaN on error. */
double
get_e(const body *b);

/** Returns the period of b. */
double
period(const body *b);

/** Constructs the coordinate system where \f$\mathbf{\hat{x}}\f$
    points toward periapse, \f$\mathbf{\hat{z}}\f$ points along the
    angular momentum, and \f$\mathbf{y} = \mathbf{z} \times
    \mathbf{x}\f$. */
void
body_coordinate_system(const body *b, double xhat[3], double yhat[3], double zhat[3]);

/** Stores the instantaneous position and velocity of the given body
    at the given eccentric anomaly in r and v.  Returns non-zero on
    error. */
int
E_to_rv(const body *b, const double E, double r[3], double v[3]);

/** Fills in the given body from the mass and orbital elements (see
 body structure for units):
 
 @param a semi-major axis.
 
 @param e eccentricity

 @param I inclination of orbital plane in degrees (I > 90 means retrograde orbit)

 @param Omega Longitude of ascending node in degrees.

 @param omega Argument of periapse in degrees.

 @param spin The spin vector *in the body frame*.  For example, spin =
 {0, 0, omega} represents the body spinning in the plane of its orbit
 (spin parallel to L).

 @param Qp is the adjusted quality factor for the tidal dissipation.

 @param inertia is the moment of inertia of the body (in units of
   m*a^2).

 @param R body radius (in units of a). */
void
init_body_from_elements(body *b,
                        const double m, const double a, const double e, const double I, 
                        const double Omega, const double omega, 
                        const double spin[3], 
                        const double Qp, const double inertia, const double R);

/** Produce a central body from the given properties. */
void
init_central_body(central_body *b, const double Qp, const double I, 
                  const double R, const double spin[3]);

/** Sets the orbital elements given a body, b. */
void
elements_from_body(const body *b,
                   double *e, double *I, double *Omega, double *omega);

/** Sets the spin to vector of the given body to synchronous rotation. */
void
body_set_synchronous_spin(body *b);

/** Computes the instantaneous gravitational acceleration on b1 due to
    b2, and the corresponding rate of change of b1's orbit. */
void
body_instantaneous_rhs(const double eps,
                       const body *b1, const double E1,
                       const body *b2, const double E2,
                       double dbdt[BODY_VECTOR_SIZE]);

#define NELEMENTS 5
#define A_INDEX 0
#define E_INDEX 1
#define I_INDEX 2
#define OMEGA_INDEX 3
#define oMEGA_INDEX 4

/** Fills deldt with [adot, edot, Idot, OmegaDot, omegaDot] from the
   given rates of change of the components of b in dbdt. */
void
body_derivs_to_orbital_elements(const body *b, const double dbdt[BODY_VECTOR_SIZE], double deldt[NELEMENTS]);

/** Computes the AMD of the given system. */
double
body_system_amd(const body bs[], const size_t n);

/* raw_average.c */

/** Numerically averaged RHS of the evolution equations for body 1 due
    to body 2. */
void
raw_average_rhs(const double eps, const body *b1, const body *b2,
                gsl_integration_workspace *ws1, const size_t ws1_limit,
                gsl_integration_workspace *ws2, const size_t ws2_limit,
                const double epsabs, const double epsrel,
                double rhs[BODY_VECTOR_SIZE]);

/* analytic_average.c */

/** Fills f with the acceleration on a body at position rp from body b
    averaged over the orbit of body b.  The softening parameter is
    eps. */
int
force_averaged_unprimed(const double eps, const double rp[3], const body *b, double f[3]);

/** Returns the time derivative of b1's components (including the mass
    and semi-major axis, which are always constant in the secular
    approximation) in rhs.  Will return GSL_FAILURE in the event of
    some failure. */
int
average_rhs(const double eps, const body *b1, const body *b2, 
            const double epsabs, double rhs[BODY_VECTOR_SIZE]);

/* advancer.c */

/** Given the number of bodies in the sysetm, returns the number of
    elements in the system vector (including the central body).  This
    is useful when allocating GSL odeiv objects, since these require a
    vector size. */
size_t
body_size_to_vector_size(const size_t nbodies);

/** Convert between arrays of bodies and arrays of doubles. */
void
bodies_to_vector(const central_body *bc, const body bs[], const size_t nbodies, double y[]);

/** Convert between arrays of doubles and arrays of bodies. */
void
vector_to_bodies(const double y[], const size_t nbodies, central_body *bc, body bs[]);

/** Advance routine.  Behaves similarly to gsl_odeiv_evolve_apply, but
   specialized to systems of bodies.  The e, con, and step arguments
   are as for gsl_odeiv_evolve_apply.  The evolve_system procedure
   uses the provided GSL objects to advance the system toward a time
   t1 (the actual time of the system after a step is returned in t,
   which is guaranteed to never pass t1).  The parameter h is the
   initial guess at a stepsize (the routine will give the actual
   stepsize chosen to maintain the accuracy requirements in the
   control object con).  bs is the initial state of the system, which
   will be overwritten as the routine progresses.  It should be of
   size nbodies.  y is a vector capable of storing the system state
   (see body_size_to_vector_size above).  y need not be initialized,
   but on output it will contain the vectorized state of bs.  The
   parameters epsabs and epsrel are integration tolerances for the
   averaging over rings.  They should be set smaller than the desired
   accuracy in the control object con.

   Typically, evolve_system is called in a loop that terminates when
   the output time is exactly t1.  In this use case, subsequent steps
   should use the parameters from the previous step (t, h, bs, y).
   When the next step is *not* a continuation of the previous step,
   remember to call the appropriate gsl_odeiv_evolve_reset procedure
   on e and gsl_odeiv_step_reset on step.  

  */
int
evolve_system(gsl_odeiv_evolve *e, gsl_odeiv_control *con, gsl_odeiv_step *step, 
              double *t, const double t1, double *h, central_body *bc, body bs[], 
              double y[], const size_t nbodies, const double epsquad, 
              const double eps);

/** A new type of control object, does the usual checks of the
    standard GSL control on the absolute error of y(t), but also
    checks the orthogonality of A and L for each body, and the
    constraint A^2 + L^2 = 1. */
gsl_odeiv_control *
gsl_odeiv_control_secular_new(double epsabs);

/* read_write.c */

/** Exects body in Runge-Lenz/L coordinates as 

   m a Qp I R L0 L1 L2 A0 A1 A2 Omega1 Omega2 Omega3

   Qp, I, and R are the dissipation constant Qp, moment of inertia,
   and radius (in units consistent with a and m).  Omega is the spin
   vector of the body, also in units consistent with a and m.
   Whitespace is ignored.  Note that these coordinates should satisfy
   the constraints A*L = 0 and A^2 + L^2 = 1, but read_body does not
   check this!  Returns number of bodies read. */
int
read_body(FILE *stream, body *b);

/** Read the central body properties from the stream in the format 

    Qp I R Omega1 Omega2 Omega3

    with the corresponding meanings for a body.
 */
int
read_central_body(FILE *stream, central_body *bc);

/** Expects elements in the following format: 

    m a e I Omega omega Qp I R S1 S2 S3

    whitespace is ignored. Qp, I, R are the dissipation constant,
    moment of inertia and radius, respectively.  S is the spin vector
    for the object.  */
int
read_body_from_elements(FILE *stream, body *b);

/** Binary version of read_body. */
int
read_body_bin(FILE *stream, body *b);

/** Binary version of read central body. */
int
read_central_body_bin(FILE *stream, central_body *bc);

/** Writes b to stream, in the format of read_body above:

   m a Qp I R L0 L1 L2 A0 A1 A2 S0 S1 S2*/
int
write_body(FILE *stream, const body *b);

/** Writes a central body to the stream in the format

   1.0 0.0 Qp I R 0.0 0.0 0.0 0.0 0.0 0.0 S0 S1 S2

 */
int 
write_central_body(FILE *stream, const central_body *b);

/** Writes the state of b to stream as 

   m a e I Omega omega Qp I R S0 S1 S2\n

   See read_body_from_elements.
*/
int
write_body_elements(FILE *stream, const body *b);

/** Writes the state of the central body to stream, in pseudo-element form:

    1.0 0.0 0.0 0.0 0.0 0.0 Qp I R S0 S1 S2
*/
int
write_central_body_elements(FILE *stream, const central_body *bc);

/** Writes b to stream in a binary format. */
int
write_body_bin(FILE *stream, const body *b);

/** Writes the central body in a binary format. */
int
write_central_body_bin(FILE *stream, const central_body *bc);

/* tides.c */

/* Computes the rate of change of the body orbital elements and the
   Sun's spin due to their tidal interaction. */
void
tidal_rhs(const body *b, const double QpSun, const double RSun, const double ISun,
          const double OmegaSun[3], double brhs[BODY_VECTOR_SIZE], double srhs[3]);
#endif /* __RINGS_H__ */
