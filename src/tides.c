/* Tidal equations taken from Barker and Ogilvie, On the Tidal
   Evolution of Hot Jupiters on Inclined Orbits, MNRAS, Vol. 395,
   Issue 4, pp 2268--2287, 2009, arXiv:0902.4563, DOI
   10.1111/j.1365-2966.2009.14694.x  */

#include"rings.h"
#include<math.h>
#include<string.h>

static double pow_int(const double x, const int n) {
  if (n < 0) {
    return pow_int(1.0/x, -n);    
  } else if (n == 0) {
    if (x == 0.0) {
      return 0.0/0.0;
    } else {
      return 1.0;
    }
  } else if (n == 1) {
    return x;
  } else if (n % 2 == 0) {
    double sqrt_xn = pow_int(x, n/2);

    return sqrt_xn*sqrt_xn;
  } else {
    double sqrt_xn = pow_int(x, (n-1)/2);

    return x * sqrt_xn * sqrt_xn;
  }
}

static double ff1(const double e2) {
  double e4 = e2*e2;
  double e6 = e4*e2;
  double sqrtDenom = sqrt(1.0 - e2);
  double denom = pow_int(sqrtDenom, 13);

  return (1.0 + 15.0/4.0*e2 + 15.0/8.0*e4 + 5.0/64.0*e6)/denom;
}

static double ff2(const double e2) {
  double e4 = e2*e2;
  double denom = pow_int(1.0 - e2, 5);

  return (1.0 + 3.0/2.0*e2 + 1.0/8.0*e4)/denom;
}

static double ff3(const double e2) {
  double e4 = e2*e2;
  double denom = pow_int(1.0 - e2, 5);

  return (1.0 + 9.0/2.0*e2 + 5.0/8.0*e4)/denom;
}

static double ff4(const double e2) {
  double e4 = e2*e2;
  double e6 = e2*e4;
  double denom = pow_int(sqrt(1.0-e2), 13);
  
  return (1.0 + 15.0/2.0*e2 + 45.0/8.0*e4 + 5.0/16.0*e6)/denom;
}

static double ff5(const double e2) {
  double denom = pow_int(1.0 - e2, 5);

  return (3.0 + 0.5*e2)/denom;
}

static double ff6(const double e2) {
  double e4 = e2*e2;
  double e6 = e4*e2;
  double e8 = e6*e2;
  double denom = pow_int(1.0-e2, 8);

  return (1.0 + 31.0/2.0*e2 + 255.0/8.0*e4 + 185.0/16.0*e6 + 25.0/64.0*e8)/denom;
}

static void
tfs(const body *b, const double QpSun, const double RSun,
    double *tb, double *tSun) {
  double n = mean_motion(b);

  *tb = 1.0 / (9.0*n/(2.0*b->Qp) * (1.0 / b->m) * pow_int(b->R / b->a, 5));
  *tSun = 1.0 / (9.0*n/(2.0*QpSun) * (b->m / 1.0) * pow_int(RSun / b->a, 5));
}

void
tidal_rhs(const body *b, const double QpSun, const double RSun, const double ISun,
          const double OmegaSun[3], double brhs[BODY_VECTOR_SIZE], double srhs[3]) {
  double e = get_e(b);
  double sqrt1me2 = sqrt(1.0 - e*e);
  double e2 = e*e;
  double n = mean_motion(b);
  double a2 = b->a*b->a;
  double mu = b->m / (1.0 + b->m);

  double Omegab2 = dot(b->spin, b->spin);
  double OmegaSun2 = dot(OmegaSun, OmegaSun);
  
  double Omegab_dot_e = dot(b->spin, b->A);
  double Omegab_dot_h = n*a2*dot(b->L, b->spin);
  
  double Omegasun_dot_e = dot(OmegaSun, b->A);
  double Omegasun_dot_h = n*a2*dot(OmegaSun, b->L);

  double h = n*a2*sqrt1me2;

  double f1 = ff1(e2);
  double f2 = ff2(e2);
  double f3 = ff3(e2);
  double f4 = ff4(e2);
  double f5 = ff5(e2);
  double f6 = ff6(e2);

  double tBody, tSun;

  double HBody, HSun;

  double adot, ndot;

  int i;

  tfs(b, QpSun, RSun, &tBody, &tSun);

  memset(brhs, 0, BODY_VECTOR_SIZE*sizeof(double));
  memset(srhs, 0, 3*sizeof(double));

  for (i = 0; i < 3; i++) {
    /* A15 from Barker and Ogilvie */
    srhs[i] = mu / (ISun*tSun) * (Omegasun_dot_e/(2.0*n)*f5*h*b->A[i] -
                                  f3*h/(2.0*n)*OmegaSun[i] + 
                                  (f4 - Omegasun_dot_h/(2.0*n*h)*f2)*n*a2*b->L[i]);
  }

  for (i = 0; i < 3; i++) {
    /* A14 from Barker and Ogilvie */
    brhs[BODY_SPIN_INDEX + i] = 
      mu / (b->I*tBody) * (Omegab_dot_e/(2.0*n)*f5*h*b->A[i] -
                           f3*h/(2.0*n)*b->spin[i] +
                           (f4 - Omegab_dot_h/(2.0*n*h)*f2)*n*a2*b->L[i]);
  }
  
  /* A27 */
  HBody = -mu*h/(n*tBody)*( 0.5*(Omegab2*f3 + Omegab_dot_h*Omegab_dot_h/(h*h)*f2 - Omegab_dot_e*Omegab_dot_e*f5) -
                            2.0*n*Omegab_dot_h*f4/h + 
                            n*n*f6 );

  /* A28 */
  HSun =  -mu*h/(n*tSun)*( 0.5*(OmegaSun2*f3 + Omegasun_dot_h*Omegasun_dot_h/(h*h)*f2 - 
                                Omegasun_dot_e*Omegasun_dot_e*f5) -
                           2.0*n*Omegasun_dot_h*f4/h +
                           n*n*f6 );

  /* From A25 and A26. */
  adot = 2.0*a2/b->m*( HBody + HSun - ISun*dot(OmegaSun, srhs) - b->I*dot(b->spin, brhs+BODY_SPIN_INDEX) );

  brhs[BODY_a_INDEX] = adot;

  ndot = -3.0/2.0*adot/b->a*n;

  for (i = 0; i < 3; i++) {
    /* A7 and A8, plus compensating for h[i] = n a^2 L[i]. */
    brhs[BODY_L_INDEX+i] = 
      (-1.0/tBody * (Omegab_dot_e/(2.0*n)*f5*h*b->A[i] - 
                     f3*h/(2.0*n)*b->spin[i] +
                     (f4 - Omegab_dot_h/(2.0*n*h)*f2)*n*a2*b->L[i]) -
       1.0/tSun  * (Omegasun_dot_e/(2.0*n)*f5*h*b->A[i] - 
                    f3*h/(2.0*n)*OmegaSun[i] +
                    (f4 - Omegasun_dot_h/(2.0*n*h)*f2)*n*a2*b->L[i]) -
       ndot*a2*b->L[i] -
       2.0*n*b->a*adot*b->L[i])/(n*a2);
  }

  for (i = 0; i < 3; i++) {
    /* A10 and A11 */
    brhs[BODY_A_INDEX+i] = 
      -1.0/tBody * (Omegab_dot_e/(2.0*n)*f2*n*a2*b->L[i] + 9.0*(f1*h - 11.0/18.0*Omegab_dot_h/n*f2)*b->A[i]) -
      1.0/tSun * (Omegasun_dot_e/(2.0*n)*f2*n*a2*b->L[i] + 9.0*(f1*h - 11.0/18.0*Omegasun_dot_h/n*f2)*b->A[i]);
    brhs[BODY_A_INDEX+i] /= h;
  }
}
