#include"rings.h"

#include<gsl/gsl_sf.h>
#include<math.h>

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
  *Bsine = a*sqrt(1.0 - e*e)*dot(rp, yhat);

  *C = a*a*e*e;
}

static void
lambda_roots(const double A, const double Bcose, const double Bsine, const double C,
             double *l0, double *l1, double *l2) {
  double CmA = C - A;
  double CmA2 = CmA*CmA;
  double CmA3 = CmA*CmA2;
  double B2 = Bsine*Bsine + Bcose*Bcose;
  double Q = 1.0/9.0*CmA2 - 1.0/3.0*(B2 - A*C);
  double R = 1.0/27.0*CmA3 - 1.0/6.0*CmA*(B2 - A*C) + 0.5*C*Bsine*Bsine;
  double sqrtQ = sqrt(Q);
  double theta = acos(R/(sqrtQ*sqrtQ*sqrtQ));

  *l0 = -2.0*sqrtQ*cos((theta + 2.0*M_PI)/3.0) - 1.0/3.0*CmA;
  *l1 = -2.0*sqrtQ*cos((theta - 2.0*M_PI)/3.0) - 1.0/3.0*CmA;
  *l2 = -2.0*sqrtQ*cos(theta/3.0) - 1.0/3.0*CmA;
}

static double
lambda_equation(const double A, const double Bcose, const double Bsine, const double C,
                const double lam) {
  double lam2 = lam*lam;
  double lam3 = lam2*lam;
  double B2 = Bcose*Bcose + Bsine*Bsine;

  return lam3 + (C - A)*lam2 + (B2 - A*C)*lam + C*Bsine*Bsine;
}

static void
Qmatrix(const double A, const double Bcose, const double Bsine, const double C,
        const double l0, const double l1, const double l2, 
        double Q[3][3]) {
  int i,j,k;

  Q[0][0] = sqrt(l0*(l0+C)/((l0-l1)*(l0-l2)));
  Q[0][1] = sqrt(l1*(l1+C)/((l0-l1)*(l1-l2)));
  Q[0][2] = sqrt(fabs(l2)*(l2+C)/((l0-l2)*(l1-l2)));
  Q[1][0] = Bsine*sqrt((l0+C)/(l0*(l0-l1)*(l0-l2)));
  Q[1][1] = Bsine*sqrt((l1+C)/(l1*(l0-l1)*(l1-l2)));
  Q[1][2] = -Bsine*sqrt((l2+C)/(fabs(l2)*(l0-l2)*(l1-l2)));
  Q[2][0] = Bcose*sqrt(l0/((l0+C)*(l0-l1)*(l0-l2)));
  Q[2][1] = Bcose*sqrt(l1/((l1+C)*(l0-l1)*(l1-l2)));
  Q[2][2] = Bcose*sqrt(fabs(l2)/((l2+C)*(l0-l2)*(l1-l2)));
}

static void
UV_from_Q(const double Q[3][3], const double e, double U[3], double V[3]) {
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

  Ek = gsl_sf_ellint_Ecomp(k, GSL_PREC_DOUBLE);
  Kk = gsl_sf_ellint_Kcomp(k, GSL_PREC_DOUBLE);

  for (i = 0; i < 3; i++) {
    f[i] = 2.0*b->m/M_PI*sqrt(l0-l2)/((l0-l1)*(l1-l2))*((k2*FU[i] + FV[i])*Ek - (1.0 - k2)*FV[i]*Kk);
  }
}

typedef struct {
  double eps;
  body *bp;
  body *b;
  int comp;
} avg_data;

static double
average_integrand(const double Ep, void *vdata) {
  avg_data *data = (avg_data *)vdata;

  if (data->comp == BODY_M_INDEX) {
    return 0.0;
  } else {
    double rp[3], vp[3];
    double f[3];
    double ep = norm(data->bp->A);
    double fac = (1.0 - ep*cos(Ep))/(2.0*M_PI);
    
    E_to_rv(data->bp, Ep, rp, vp);
    force_averaged_unprimed(data->eps, rp, data->b, f);
    
    if (data->comp == BODY_a_INDEX) {
      double np = mean_motion(data->bp);
      double ap = data->bp->a;

      return fac*2.0*dot(vp, f)/(np*np*ap);
    } else if (BODY_A_INDEX <= data->comp && data->comp < BODY_A_INDEX+3) {
      double fdv = dot(f, vp);
      double rdf = dot(rp, f);
      double rdv = dot(rp, vp);
      int i = data->comp - BODY_A_INDEX;

      return fac/(1.0 + data->bp->m)*(2.0*fdv*rp[i] - rdf*vp[i] - rdv*f[i]);
    } else if (BODY_L_INDEX <= data->comp && data->comp < BODY_L_INDEX+3) {
      double np = mean_motion(data->bp);
      double ap = data->bp->a;
      double rxf[3];
      int i = data->comp - BODY_L_INDEX;

      cross(rp, f, rxf);

      return fac*rxf[i]/(np*ap*ap);
    } else {
      fprintf(stderr, "Unknown index for RHS component in average_integrand.\n");
      abort();
    }
  }  
}

void
average_rhs(const double eps, const body *b1, const body *b2, 
            gsl_integration_workspace *ws, size_t ws_size,
            const double epsabs, const double epsrel, double rhs[BODY_VECTOR_SIZE]) {
  avg_data data;
  int i;  

  data.b = b2;
  data.bp = b1;
  data.eps = eps;

  for (i = 0; i < BODY_VECTOR_SIZE; i++) {
    gsl_function f;
    double result, err;

    data.comp = i;

    f.function = (double (*)(const double, void *))average_integrand;
    f.params = &data;

    gsl_integration_qag(&f, 0.0, 2.0*M_PI, epsabs, epsrel, ws_size, GSL_INTEG_GAUSS61, ws, &result, &err);

    rhs[i] = result;
  }
}
