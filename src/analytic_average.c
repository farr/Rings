#include"rings.h"

#include<gsl/gsl_sf.h>
#include<gsl/gsl_poly.h>
#include<gsl/gsl_integration.h>

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

  double theta = acos(cosTheta);
  
  *l0 = -2.0*sqrtQ*cos(theta/3.0 + 2.0*M_PI/3.0) - CmAO3;
  *l1 = -2.0*sqrtQ*cos(theta/3.0 - 2.0*M_PI/3.0) - CmAO3;
  *l2 = -2.0*sqrtQ*cos(theta/3.0) - CmAO3;

  if (C < 0.0) {
    fprintf(stderr, "ERROR: bad C value in lambda_roots.\n");
    exit(1);
  } else if (*l2 < -C) {
    /* Fix up assuming that the other roots are correct. */
    *l2 = Bcose*Bcose*C/((*l0+C)*(*l1+C)) - C;
  } else if (*l2 > 0.0) {
    *l2 = -Bsine*Bsine*C / ((*l0)*(*l1));
  } else if (*l1 < 0.0) {
    fprintf(stderr, "ERROR: l1 < 0 in lambda_roots.\n");
    exit(1);
  } else if (*l1 > *l0) {
    fprintf(stderr, "ERROR: l1 > l0 in lambda_roots.\n");
    exit(1);
  } else if (A < 0.0) {
    fprintf(stderr, "ERROR: A < 0 in lambda_roots.\n");
    exit(1);
  } else if (*l0 > A) {
    fprintf(stderr, "ERROR: l0 > A in lambda_roots.\n");
    exit(1);
  }
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
        const double e, 
        double Q[3][3]) {
  double l1 = fabs(l1_in);  /* l1 should always be nonnegative, but
                               roundoff.... */
  double l2 = -fabs(l2_in);  /* l2 should always be nonpositive, but
                                roundoff.... */
  
  Q[0][0] = sqrt(l0*(l0+C)/((l0-l1)*(l0-l2)));
  Q[0][1] = sqrt(l1*(l1+C)/((l0-l1)*(l1-l2)));
  Q[0][2] = sqrt(fabs(l2)*fabs(l2+C)/((l0-l2)*(l1-l2)));
  Q[1][0] = Bsine*sqrt((l0+C)/(l0*(l0-l1)*(l0-l2)));
  Q[1][1] = Bsine*sqrt((l1+C)/(l1*(l0-l1)*(l1-l2)));
  Q[1][2] = -Bsine*sqrt(fabs(l2+C)/(fabs(l2)*(l0-l2)*(l1-l2)));
  Q[2][0] = Bcose*sqrt(l0/((l0+C)*(l0-l1)*(l0-l2)));
  Q[2][1] = Bcose*sqrt(l1/((l1+C)*(l0-l1)*(l1-l2)));
  Q[2][2] = Bcose*sqrt(fabs(l2)/(fabs(l2+C)*(l0-l2)*(l1-l2)));

  if (e*e < 1e-8) {
    /* Small e => small C => small l2, correct Q[1][2] and Q[2][2] */
    Q[1][2] = - sign(Bsine)*fabs(Bcose)/sqrt(fabs((l0-l2)*(l1-l2)));
    Q[2][2] = sign(Bcose)*fabs(Bsine)/sqrt(fabs((l0-l2)*(l1-l2)));
  } else if (fabs(1.0 - e*e) < 1e-8) {
    /* Large e => small B*sin(eps) => small l1, correct Q[1][1]. */
    Q[1][1] = 1.0 - (C*l0 + C*l2 + l0*l2)*Bsine*Bsine/(2.0*l0*l0*l2*l2);
  } else if (fabs(2.0*l2/(l1 + l0)) < 1e-8) {
    /* l2 driven small => BSine similarly small, but, unlike e --> 0,
       l2 not close to C. */
    Q[1][2] = -1.0 + Bsine*Bsine*(l0*l1 + C*(l0+l1))/(2.0*l0*l0*l1*l1);
  } else if (fabs(2.0*l1/(l0+l2)) < 1e-8) {
    /* l1 driven small => BSine similarly small.*/
    Q[1][1] = 1.0 - (C*l0 + C*l2 + l0*l2)*Bsine*Bsine/(2.0*l0*l0*l2*l2);
  } else if (fabs(2.0*(l2+C)/(l1+C + l0+C)) < 1e-8) {
    /* l2 driven toward -C => BCose similarly small. */
    Q[2][2] = sqrt(fabs(l2*(l0+C)*(l1+C)/(C*(l0-l2)*(l1-l2))));
  }

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
UV_from_lambda(const double l0, const double l1, const double l2, const double A, 
               const double Bsine, const double Bcose, const double C, const double e,
               double U[3], double V[3]) {
  U[0] = ((l0*(C-Bcose*e+l0))/(l0-l1) + (l2*(C-Bcose*e+l2))/(l2-l1))/(l0-l2);
  U[1] = (Bsine*(l0*l1 - 2.0*l0*l2 + l1*l2 - C*(l0 - 2*l1 + l2) + Bcose*e*(l0 - 2*l1 + l2)))/((l0-l1)*(l0-l2)*(l1-l2));
  U[2] = (Bcose*((l0*(C-Bcose*e+l0))/((l0+C)*(l0-l1)) + (fabs(l2)*(C-Bcose*e+l2))/((l1-l2)*(l2+C))))/(l0-l2);
  V[0] = (l1*(C-Bcose*e+l1)*(l0-l2) + (l0-l1)*l2*(C-Bcose*e+l2))/((l0-l1)*(l0-l2)*(l1-l2));
  V[1] = (Bsine*(l0*l1 + C*(2.0*l0-l1-l2)+l0*l2-2.0*l1*l2+Bcose*e*(l1+l2-2.0*l0)))/((l0-l1)*(l0-l2)*(l1-l2));
  V[2] = Bcose*(l1*(C-Bcose*e+l1)/((l0-l1)*(C+l1)*(l1-l2)) + l2*(C-Bcose*e+l2)/((C+l2)*(l2-l0)*(l2-l1)));
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
    F1[i] = a*sqrt(fabs(1.0-e*e))*yhat[i];
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
  double e = norm(b->A);
  double k, k2;
  double Ek, Kk;
  int i;

  get_ABC(eps, rp, b, &A, &Bcose, &Bsine, &C);
  lambda_roots(A, Bcose, Bsine, C, &l0, &l1, &l2);
  /* Qmatrix(A, Bcose, Bsine, C, l0, l1, l2, e, Q);

     UV_from_Q(Q, norm(b->A), U, V); */

  UV_from_lambda(l0, l1, l2, A, Bsine, Bcose, C, e, U, V);
  
  get_F(rp, b, F0, F1, F2);

  get_FU_FV(U, V, F0, F1, F2, FU, FV);

  k2 = (l1-l2)/(l0-l2);

  k = sqrt(fabs(k2));

  if (isnan(k)) {
    fprintf(stderr, "k = NaN in %s, line %d\n", __FILE__, __LINE__);
    exit(1);
  }

  Ek = gsl_sf_ellint_Ecomp(k, GSL_PREC_DOUBLE);
  Kk = gsl_sf_ellint_Kcomp(k, GSL_PREC_DOUBLE);

  for (i = 0; i < 3; i++) {
    f[i] = 2.0*b->m/M_PI*sqrt(fabs(l0-l2))/((l0-l1)*(l1-l2))*((k2*FU[i] + FV[i])*Ek - (1.0 - k2)*FV[i]*Kk);
    assert(!isnan(f[i]));
  }
}

static void 
force_to_components(const double r[3], const double v[3], const double f[3], 
                    double *R, double *T, double *W) {
  double l[3];
  double rhat[3], that[3], lhat[3];

  cross(r, v, l);
  unitize(l, lhat);

  unitize(r, rhat);

  cross(lhat, rhat, that);

  *R = dot(rhat, f);
  *T = dot(that, f);
  *W = dot(lhat, f);
}

static void
force_components_integrand(const double eps, const body *b, const body *bp, const double E,
                           double *R, double *T, double *W) {
  double r[3], v[3], f[3];

  E_to_rv(bp, E, r, v);
  force_averaged_unprimed(eps, r, b, f);
  force_to_components(r, v, f, R, T, W);
}

static void 
average_force_components(const double eps, const body *b, const body *bp, 
                         const double epsacc,
                         double Rc[3], double Rs[2],
                         double Tc[3], double Ts[2], 
                         double Wc[3], double Ws[2]) {
  double rp[3], vp[3];
  double ep = norm(bp->A);
  double a = bp->a, n = mean_motion(bp);
  double R, T, W;
  double err;
  int N, i;
  const int NMax = 1048576;

  memset(Rc, 0, 3*sizeof(double));
  memset(Rs, 0, 2*sizeof(double));
  memset(Tc, 0, 3*sizeof(double));
  memset(Ts, 0, 2*sizeof(double));
  memset(Wc, 0, 3*sizeof(double));
  memset(Ws, 0, 2*sizeof(double));

  /* Compute the E = 0 term. */
  force_components_integrand(eps, b, bp, 0.0, &R, &T, &W);
  for (i = 0; i < 3; i++) {
    Rc[i] += R;
    Tc[i] += T;
    Wc[i] += W;
  }

  N = 1;
  do {
    double h;

    N *= 2;
    h = 1.0 / N;

    for (i = 0; i < 3; i++) {
      Rc[i] /= 2.0;
      Tc[i] /= 2.0;
      Wc[i] /= 2.0;
    }

    for (i = 0; i < 2; i++) {
      Rs[i] /= 2.0;
      Ts[i] /= 2.0;
      Ws[i] /= 2.0;
    }

    for (i = 1; i < N; i += 2) {
      double E = 2.0*M_PI*((double) i)/((double) N);
      int j;

      force_components_integrand(eps, b, bp, E, &R, &T, &W);

      Rc[0] += h*R;
      Tc[0] += h*T;
      Wc[0] += h*W;
      for (j = 0; j < 2; j++) {
        double c = cos((j+1)*E);
        double s = sin((j+1)*E);

        Rc[j+1] += h*R*c;
        Tc[j+1] += h*T*c;
        Wc[j+1] += h*W*c;

        Rs[j] += h*R*s;
        Ts[j] += h*T*s;
        Ws[j] += h*W*s;
      }
    }

    err = fabs(ep*Rs[0] + sqrt(1.0 - ep*ep)*Tc[0])/(a*n*n);
  } while ((N < 16 || err > epsacc) && (N < NMax));

  if (N >= NMax) {
    int i;
    fprintf(stderr, "ERROR: max iters exceeded in numerical averaging (file %s, line %d, %s)\n",
            __FILE__, __LINE__, __func__);
    fprintf(stderr, "ERROR: err = %g\n", err);
    fprintf(stderr, "ERROR: Rc[0] = %g, Tc[0] = %g, Wc[0] = %g\n",
            Rc[0], Tc[0], Wc[0]);
    for (i = 0; i < 2; i++) {
      fprintf(stderr, "ERROR: Rc[%d] = %g, Tc[%d] = %g, Wc[%d] = %g\n",
              i+1, Rc[i+1], i+1, Tc[i+1], i+1, Wc[i+1]);
      fprintf(stderr, "ERROR: Rs[%d] = %g, Ts[%d] = %g, Ws[%d] = %g\n",
              i+1, Rs[i], i+1, Ts[i], i+1, Ws[i]);
    }
    exit(1);
  }
}

void
average_rhs(const double eps, const body *b1, const body *b2, 
            const double epsabs, double rhs[BODY_VECTOR_SIZE]) {
  int i;
  double Rc[3], Tc[3], Wc[3];
  double Rs[2], Ts[2], Ws[2];
  double a = b1->a, e = norm(b1->A);
  double e2 = e*e;
  double sqrt1me2 = norm(b1->L);
  double N[3], A[3];
  double xhat[3], yhat[3], zhat[3];
  double n = mean_motion(b1);
  double na2 = n*a*a;
  double na = n*a;

  average_force_components(eps, b2, b1, epsabs, Rc, Rs, Tc, Ts, Wc, Ws);

  body_coordinate_system(b1, xhat, yhat, zhat);

  N[0] = a*sqrt1me2*(Ws[0] - 0.5*e*Ws[1]);
  N[1] = -a*((1.0+e2)*Wc[1] - 1.5*e*Wc[0] - 0.5*e*Wc[2]);
  N[2] = a*((1.0+0.5*e2)*Tc[0] - 2.0*e*Tc[1] + 0.5*e2*Tc[2]);

  A[0] = sqrt1me2/(2.0*na)*((4.0*Tc[1] - e*Tc[2] - 3.0*e*Tc[0]) + 2.0*sqrt1me2*Rs[0]);
  A[1] = 1.0/(2.0*na)*(2.0*(2.0-e2)*Ts[0] - e*Ts[1] - 2.0*sqrt1me2*(Rc[1] - e*Rc[0]));
  A[2] = -e/na*(Ws[0] - 0.5*e*Ws[1]);

  memset(rhs, 0, BODY_VECTOR_SIZE*sizeof(double));

  for (i = 0; i < 3; i++) {
    rhs[BODY_L_INDEX + i] = (N[0]*xhat[i] + N[1]*yhat[i] + N[2]*zhat[i])/na2;
    rhs[BODY_A_INDEX + i] =  A[0]*xhat[i] + A[1]*yhat[i] + A[2]*zhat[i];
  }
}
