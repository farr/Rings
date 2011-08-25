#include"rings.h"

#include<gsl/gsl_sf.h>
#include<gsl/gsl_poly.h>

#include<math.h>
#include<assert.h>
#include<string.h>

static void
get_ABC(const double eps, const double rp[3], const body *b, 
        double *A, double *Bcose, double *Bsine, double *C) {
  double rpmag = norm(rp);
  double e = norm(b->A);
  double xhat[3], yhat[3], zhat[3];
  double a = b->a;

  body_coordinate_system(b, xhat, yhat, zhat);

  *A = rpmag*rpmag + a*a + eps*eps + 2.0*a*e*dot(rp, xhat);

  *Bcose = a*dot(rp, xhat) + a*a*e;
  *Bsine = a*sqrt(fabs(1.0 - e*e))*dot(rp, yhat);

  *C = a*a*e*e;
}

static void
lambda_roots(const double A, const double Bcose, const double Bsine, const double C,
             double *l0, double *l1, double *l2) {
  int nroots;
  double B2 = Bcose*Bcose + Bsine*Bsine;
  double CmA = C - A;
  double CmA2 = CmA*CmA;
  double CmA3 = CmA2*CmA;
  double disc = B2 - A*C;
  double Q = 1.0/9.0*CmA2 - 1.0/3.0*disc;
  double R = 1.0/27.0*CmA3 - 1.0/6.0*CmA*disc + 0.5*Bsine*Bsine*C;
  double sqrtQ = sqrt(Q);
  double CmAO3 = CmA/3.0;

  double cosTheta = R / (sqrtQ*sqrtQ*sqrtQ);

  if (cosTheta > 1.0) cosTheta = 1.0;
  if (cosTheta < -1.0) cosTheta = -1.0;

  double theta = acos(cosTheta);
  
  *l0 = -2.0*sqrtQ*cos(theta/3.0 + 2.0*M_PI/3.0) - CmAO3;
  *l1 = -2.0*sqrtQ*cos(theta/3.0 - 2.0*M_PI/3.0) - CmAO3;
  *l2 = -2.0*sqrtQ*cos(theta/3.0) - CmAO3;
}

static double
lambda_equation(const double A, const double Bcose, const double Bsine, const double C,
                const double lam) {
  double lam2 = lam*lam;
  double lam3 = lam2*lam;
  double B2 = Bcose*Bcose + Bsine*Bsine;

  return lam3 + (C - A)*lam2 + (B2 - A*C)*lam + C*Bsine*Bsine;
}

static double
sign(const double x) { if (x < 0.0) { return -1.0; } else { return 1.0; } }

static void
Qmatrix(const double A, const double Bcose, const double Bsine, const double C,
        const double l0, const double l1_in, const double l2_in, 
        double Q[3][3]) {
  int i,j,k;
  double l1 = fabs(l1_in);  /* l1 should always be nonnegative, but
                               roundoff.... */
  double l2 = -fabs(l2_in);  /* l2 should always be nonpositive, but
                                roundoff.... */

  double root_sum = l0+l1+fabs(l2);
  double small = sqrt(DBL_EPSILON)*root_sum;

  /* We have several issues to worry about here, hence the careful
     finessing with fabs(...) and re-writing of some terms relative to
     T&T. 

     As e --> 0, C --> 0, and -C < l2 < 0 means that l2 becomes
     "squeezed" between -C and 0.  In this case l2/(l2+C) --> 0/0.
     So, we want to avoid this ratio.
       
     As e --> 1, 0 < l1 < a^2(1-e^2) gets squeezed toward zero.  In
     this circumstance, the l1 in the denominator Q[1][1] is
     problematic.  In fact, as e --> 1, Bsine --> 0 as well, so we
     again have a 0/0 situation.

     When rp*yhat --> 0, Bsine --> 0 (this happens at least once an
     orbit); because l0*l1*l2 = -B^2 C Sin^2 epsilon, at least one of
     l2 or l1 is driven to zero.  

     We can fix these problems by using the relations

     l0*l1*l2 = -Bsine^2 C
     (l0+C)*(l1+C)*(l2+C) = Bcose^2 C

     which imply that 

     (l2+C)/l2 = cose^2/sine^2 l0*l1/((l0+C)*(l1+C))

     to remove l2, l1, and (l2+C) terms from the denominator of the Q
     components.  The result appears below.  We also added fabs(...)
     to l1 and (l2+C) because they should be positive, but may be
     negative due to roundoff.  (fabs(l2) already appears because l2
     *should* be negative.)

  */

  /* Though l2 + C should always be positive, sometimes it gets small
     enough that roundoff matters; in this case, we need the
     fabs(l2+C) calls. */
  Q[0][0] = sqrt(l0*(l0+C)/((l0-l1)*(l0-l2)));
  Q[0][1] = sqrt(l1*(l1+C)/((l0-l1)*(l1-l2)));
  Q[0][2] = sqrt(fabs(l2)*fabs(l2+C)/((l0-l2)*(l1-l2)));
  Q[1][0] = Bsine*sqrt((l0+C)/(l0*(l0-l1)*(l0-l2)));
  Q[1][1] = Bsine*sqrt((l1+C)/(l1*(l0-l1)*(l1-l2)));
  Q[1][2] = -Bsine*sqrt(fabs(l2+C)/(fabs(l2)*(l0-l2)*(l1-l2)));
  Q[2][0] = Bcose*sqrt(l0/((l0+C)*(l0-l1)*(l0-l2)));
  Q[2][1] = Bcose*sqrt(l1/((l1+C)*(l0-l1)*(l1-l2)));
  Q[2][2] = Bcose*sqrt(fabs(l2)/(fabs(l2+C)*(l0-l2)*(l1-l2)));
}

static void
UV_from_Q(double Q[3][3], const double e, double U[3], double V[3]) {
  U[0] = Q[0][0]*Q[0][0] - e*Q[0][0]*Q[2][0] + Q[0][2]*Q[0][2] - e*Q[0][2]*Q[2][2];
  U[1] = Q[0][0]*Q[1][0] - e*Q[1][0]*Q[2][0] + Q[0][2]*Q[1][2] - e*Q[1][2]*Q[2][2];
  U[2] = Q[0][0]*Q[2][0] - e*Q[2][0]*Q[2][0] + Q[0][2]*Q[2][2] - e*Q[2][2]*Q[2][2];
  
  V[0] = Q[0][1]*Q[0][1] - e*Q[0][1]*Q[2][1] - Q[0][2]*Q[0][2] + e*Q[0][2]*Q[2][2];
  V[1] = Q[0][1]*Q[1][1] - e*Q[1][1]*Q[2][1] - Q[0][2]*Q[1][2] + e*Q[1][2]*Q[2][2];
  V[2] = Q[0][1]*Q[2][1] - e*Q[2][1]*Q[2][1] - Q[0][2]*Q[2][2] + e*Q[2][2]*Q[2][2];
}

static void
get_F(const double rp[3], const body *b, double F0[3], double F1[3], double F2[3]) {
  double xhat[3], yhat[3], zhat[3];
  double a = b->a;
  double e = norm(b->A);
  int i;

  body_coordinate_system(b, xhat, yhat, zhat);

  for (i = 0; i < 3; i++) {
    F0[i] = -rp[i] - a*e*xhat[i];
    F1[i] = a*sqrt(1.0-e*e)*yhat[i];
    F2[i] = a*xhat[i];
  }
}

static void
get_FU_FV(const double U[3], const double V[3], 
          const double F0[3], const double F1[3], const double F2[3],
          double FU[3], double FV[3]) {
  int i;

  for (i = 0; i < 3; i++) {
    int j;

    FU[i] = U[0]*F0[i] + U[1]*F1[i] + U[2]*F2[i];
    FV[i] = V[0]*F0[i] + V[1]*F1[i] + V[2]*F2[i];
  }
}

void
force_averaged_unprimed(const double eps, const double rp[3], const body *b, double f[3]) {
  double l0, l1, l2;
  double Q[3][3];
  double U[3], V[3];
  double F0[3], F1[3], F2[3];
  double FU[3], FV[3];
  double A, Bsine, Bcose, C;
  double k, k2;
  double Ek, Kk;
  int i;

  get_ABC(eps, rp, b, &A, &Bcose, &Bsine, &C);
  lambda_roots(A, Bcose, Bsine, C, &l0, &l1, &l2);
  Qmatrix(A, Bcose, Bsine, C, l0, l1, l2, Q);

  UV_from_Q(Q, norm(b->A), U, V);
  
  get_F(rp, b, F0, F1, F2);

  get_FU_FV(U, V, F0, F1, F2, FU, FV);

  k2 = (l1-l2)/(l0-l2);
  k = sqrt(fabs(k2));

  assert(!(isnan(k))); /* Workaround for infinite loop bug in GSL ellint. */

  Ek = gsl_sf_ellint_Ecomp(k, GSL_PREC_DOUBLE);
  Kk = gsl_sf_ellint_Kcomp(k, GSL_PREC_DOUBLE);

  for (i = 0; i < 3; i++) {
    f[i] = 2.0*b->m/M_PI*sqrt(l0-l2)/((l0-l1)*(l1-l2))*((k2*FU[i] + FV[i])*Ek - (1.0 - k2)*FV[i]*Kk);
  }
}

typedef struct {
  double eps;
  const body *bp;
  const body *b;
} avg_data;

static int
average_integrand(const double Ep, void *vdata, double result[BODY_VECTOR_SIZE]) {
  avg_data *data = (avg_data *)vdata;
  double rp[3], vp[3];
  double f[3];
  double ep = norm(data->bp->A);
  double fac = (1.0 - ep*cos(Ep))/(2.0*M_PI);
  double fdv;
  double rdf;
  double rdv;
  double rxf[3];
  double ap = data->bp->a;
  double np = mean_motion(data->bp);
  int i;

  result[BODY_M_INDEX] = 0.0;
  result[BODY_Qp_INDEX] = 0.0;
  result[BODY_I_INDEX] = 0.0;
  result[BODY_R_INDEX] = 0.0;

  E_to_rv(data->bp, Ep, rp, vp);

  force_averaged_unprimed(data->eps, rp, data->b, f);
  
  fdv = dot(f, vp);
  rdf = dot(rp, f);
  rdv = dot(rp, vp);
  cross(rp, f, rxf);

  result[BODY_a_INDEX] = 2.0*fac/(np*np*ap)*fdv;
  
  for (i = 0; i < 3; i++) {
    result[BODY_A_INDEX+i] = fac/(1.0 + data->bp->m)*(2.0*fdv*rp[i] - rdf*vp[i] - rdv*f[i]);
  }

  for (i = 0; i < 3; i++) {
    result[BODY_L_INDEX+i] = fac*rxf[i]/(np*ap*ap);
  }

  memset(result+BODY_SPIN_INDEX, 0, 3*sizeof(double));

  return GSL_SUCCESS;
}

int
average_rhs(const double eps, const body *b1, const body *b2, 
            const double epsabs, double rhs[BODY_VECTOR_SIZE]) {
  avg_data data;

  data.b = b2;
  data.bp = b1;
  data.eps = eps;

  return quad(average_integrand, &data, 0.0, 2.0*M_PI, epsabs, rhs);
}
