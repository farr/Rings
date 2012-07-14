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

static void 
compute_XYZVW(const double m1, const double tV, const double k, const double Roa, const double n, const double a,  
              const double e, const double m2, const double xhat[3], const double yhat[3], 
              const double zhat[3], const double spin[3],
              double *X, double *Y, double *Z, double *V, double *W) {
  double e2 = e*e;
  double e4 = e2*e2;
  double e6 = e4*e2;

  double Roa5 = pow_int(Roa, 5);
  double aoR8 = pow_int(Roa, -8);
  
  double mu = m1*m2/(m1+m2);
  
  double tF = tV/9.0*aoR8*m1*m1/(m1+m2)/m2*pow_int(1.0 + 2.0*k, -2);

  double Ox = dot(spin, xhat);
  double Oy = dot(spin, yhat);
  double Oz = dot(spin, zhat);

  *X = -m2*k*Roa5/(mu*n)*Oz*Ox/pow_int(1-e2, 2) - Oy/(2*n*tF)*(1.0 + (9.0/2.0)*e2 + (5.0/8.0)*e4)/pow_int(1-e2, 5);
  *Y = -m2*k*Roa5/(mu*n)*Oz*Oy/pow_int(1-e2, 2) + Ox/(2*n*tF)*(1.0 + (3.0/2.0)*e2 + (1.0/8.0)*e4)/pow_int(1-e2, 5);
  *Z =  m2*k*Roa5/(mu*n)*((2*Oz*Oz - Oy*Oy - Ox*Ox)/(2.0*pow_int(1-e2,2)) + 15.0*m2/pow_int(a,3)*(1.0 + (3.0/2.0)*e2 + (1.0/8.0)*e4)/pow_int(1-e2,5));

  *V = 9.0/tF*((1.0 + (15.0/4.0)*e2 + (15.0/8.0)*e4 + (5.0/64.0)*e6)/sqrt(pow_int(1-e2,13)) - 11.0*Oz/(18.0*n)*(1.0 + (3.0/2.0)*e2 + (1.0/8.0)*e4)/pow_int(1-e2, 5));
  *W = 1.0/tF*((1.0 + (15.0/2.0)*e2 + (45.0/8.0)*e4 + (5.0/16.0)*e6)/sqrt(pow_int(1-e2,13)) - Oz/n*(1.0 + 3.0*e2 + (3.0/8.0)*e4)/pow_int(1-e2, 5));
}

/* Equations from Fabrycky & Tremaine (2007). */
void
tidal_rhs(const body *b, const central_body *bc, double brhs[BODY_VECTOR_SIZE], double srhs[3]) {
  double xhat[3], yhat[3], zhat[3];
  double e, n;
  
  double X1, Y1, Z1, V1, W1;
  double X2, Y2, Z2, V2, W2;
  
  double h, adot, mu, srtma, a2, a;

  int i;

  if (b->R <= 0.0 && bc->R <= 0.0) {
    /* Fast path: if R == 0, no tides! */
    memset(brhs, 0, BODY_VECTOR_SIZE*sizeof(double));
    memset(srhs, 0, 3*sizeof(double));

    return;
  }

  body_coordinate_system(b, xhat, yhat, zhat);

  e = get_e(b);
  n = mean_motion(b);
  a = b->a;
  a2 = a*a;

  h = norm(b->L);

  mu = b->m/(1.0+b->m);

  srtma = sqrt((1.0+b->m)*b->a);

  compute_XYZVW(1.0, bc->tV, bc->k, bc->R/b->a, n, b->a, e, b->m, xhat, yhat, zhat, bc->spin, &X1, &Y1, &Z1, &V1, &W1);
  compute_XYZVW(b->m, b->tV, b->k, b->R/b->a, n, b->a, e, 1.0, xhat, yhat, zhat, b->spin, &X2, &Y2, &Z2, &V2, &W2);

  memset(brhs, 0, BODY_VECTOR_SIZE*sizeof(double));
  memset(srhs, 0, 3*sizeof(double));

  adot = -2*b->a*(W1 + W2 + e*e/(1-e*e)*(V1 + V2));
  brhs[BODY_a_INDEX] = adot;

  for (i = 0; i < 3; i++) {
    brhs[BODY_A_INDEX + i] = e*((Z1 + Z2)*yhat[i] - (Y1 + Y2)*zhat[i] - (V1 + V2)*xhat[i]);
    brhs[BODY_L_INDEX + i] = h*((Y1 + Y2)*xhat[i] - (X1 + X2)*yhat[i] - (W1 + W2)*zhat[i]) - 0.5*(adot/b->a)*b->L[i];
    brhs[BODY_SPIN_INDEX + i] = mu*h*srtma/b->I*(-Y1*xhat[i] + X1*yhat[i] + W1*zhat[i]);
    srhs[i] = mu*h*srtma/bc->I*(-Y2*xhat[i] + X2*yhat[i] + W2*zhat[i]);
  }
}
